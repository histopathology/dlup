# Copyright (c) dlup contributors
"""
Classes to write image and mask files
"""
from __future__ import annotations

import abc
import pathlib
import shutil
import tempfile
from enum import Enum
from typing import Any, Generator, Iterator

import numpy as np
import numpy.typing as npt
import PIL.Image
import PIL.ImageColor
from tifffile import tifffile

import dlup
from dlup._libtiff_tiff_writer import LibtiffTiffWriter
from dlup._types import PathLike
from dlup.tiling import Grid, GridOrder, TilingMode
from dlup.utils.tifffile_utils import get_tile


class TiffCompression(str, Enum):
    """Compression types for tiff files."""

    NONE = "NONE"  # No compression
    CCITTFAX4 = "CCITTFAX4"  # Fax4 compression
    JPEG = "JPEG"  # Jpeg compression
    DEFLATE = "DEFLATE"  # zip compression
    PACKBITS = "PACKBITS"  # packbits compression
    LZW = "LZW"  # LZW compression, not implemented in tifffile
    WEBP = "WEBP"  # WEBP compression
    ZSTD = "ZSTD"  # ZSTD compression
    JP2K = "JP2K"  # JP2K compression
    JP2K_LOSSY = "JP2K_LOSSY"
    PNG = "PNG"


# Mapping to map TiffCompression to their respective values in tifffile.
TIFFFILE_COMPRESSION = {
    "NONE": None,
    "CCITTFAX4": "CCITT_T4",
    "JPEG": "jpeg",
    "DEFLATE": "deflate",
    "PACKBITS": "packbits",
    "LZW": "lzw",
    "WEBP": "webp",
    "ZSTD": "zstd",
    "JP2K": "jpeg2000",
    "JP2K_LOSSY": "jpeg_2000_lossy",
    "PNG": "png",
}


def _color_dict_to_color_lut(color_map: dict[int, str]) -> npt.NDArray[np.uint16]:
    """
    Convert a color map to a color look-up table (LUT).

    Parameters
    ----------
    color_map : dict
        Color map to convert. The keys are the indices and the values are the color names.

    Returns
    -------
    npt.NDArray[np.uint16]
        Color LUT as a 3x256 array.
    """
    # Convert color names to RGB values
    rgb_color_map = {index: PIL.ImageColor.getrgb(color_name) for index, color_name in color_map.items()}

    # Initialize a 3x256 CLUT (for 8-bit images)
    color_lut = np.zeros((3, 256), dtype=np.uint16)

    # Prepare indices and corresponding colors for assignment
    indices = np.array(list(rgb_color_map.keys()))
    colors = np.array(list(rgb_color_map.values())) * 256  # Scale to 16-bit color depth

    # Assign colors to clut using advanced indexing
    color_lut[0, indices] = colors[:, 0]  # Red channel
    color_lut[1, indices] = colors[:, 1]  # Green channel
    color_lut[2, indices] = colors[:, 2]  # Blue channel

    return color_lut


class ImageWriter(abc.ABC):
    """Base writer class"""

    def __init__(
        self,
        filename: PathLike,
        size: tuple[int, int] | tuple[int, int, int],
        mpp: float | tuple[float, float],
        tile_size: tuple[int, int] = (512, 512),
        pyramid: bool = False,
        colormap: dict[int, str] | None = None,
        compression: TiffCompression | None = TiffCompression.JPEG,
        is_mask: bool = False,
        quality: int | None = 100,
        metadata: dict[str, str] | None = None,
    ):

        if compression is None:
            compression = TiffCompression.NONE

        self._filename = filename
        self._tile_size = tile_size
        self._size = (*size[::-1], 1) if len(size) == 2 else (size[1], size[0], size[2])
        self._mpp: tuple[float, float] = (mpp, mpp) if isinstance(mpp, (int, float)) else mpp
        self._pyramid = pyramid
        self._colormap = _color_dict_to_color_lut(colormap) if colormap is not None else None
        self._compression = compression
        self._is_mask = is_mask
        self._quality = quality
        self._metadata = metadata

    def from_pil(self, pil_image: PIL.Image.Image) -> None:
        """
        Create tiff image from a PIL image

        Parameters
        ----------
        pil_image : PIL.Image
        """
        if not np.all(np.asarray(pil_image.size)[::-1] >= self._tile_size):
            raise RuntimeError(
                f"PIL Image must be larger than set tile size. Got {pil_image.size} and {self._tile_size}."
            )
        iterator = _tiles_iterator_from_pil_image(pil_image, self._tile_size, order="F")
        self.from_tiles_iterator(iterator)

    @abc.abstractmethod
    def from_tiles_iterator(self, iterator: Iterator[npt.NDArray[np.int_]]) -> None:
        """"""


class LibtiffImageWriter(ImageWriter):
    """Image writer that writes tile-by-tile to tiff using LibtiffWriter."""

    def __init__(
        self,
        filename: PathLike,
        size: tuple[int, int] | tuple[int, int, int],
        mpp: float | tuple[float, float],
        tile_size: tuple[int, int] = (512, 512),
        pyramid: bool = False,
        colormap: dict[int, str] | None = None,
        compression: TiffCompression | None = TiffCompression.JPEG,
        is_mask: bool = False,
        quality: int | None = 100,
        metadata: dict[str, str] | None = None,
    ):
        super().__init__(
            filename,
            size,
            mpp,
            tile_size,
            pyramid,
            colormap,
            compression,
            is_mask,
            quality,
            metadata,
        )

        compression_value: str
        if isinstance(self._compression, TiffCompression):
            compression_value = self._compression.value
        else:
            compression_value = self._compression

        self._writer = LibtiffTiffWriter(
            self._filename,
            self._size,
            self._mpp,
            self._tile_size,
            compression_value,
            self._quality if self._quality is not None else 100,
        )

    def from_tiles_iterator(self, iterator: Iterator[npt.NDArray[np.int_]]) -> None:
        tiles_per_row = (self._size[1] + self._tile_size[1] - 1) // self._tile_size[1]

        for idx, tile in enumerate(iterator):
            row = (idx // tiles_per_row) * self._tile_size[0]
            col = (idx % tiles_per_row) * self._tile_size[1]
            self._writer.write_tile(tile, row, col)

        if self._pyramid:
            self._writer.write_pyramid()

        self._writer.finalize()


class TifffileImageWriter(ImageWriter):
    """Image writer that writes tile-by-tile to tiff."""

    def __init__(
        self,
        filename: PathLike,
        size: tuple[int, int] | tuple[int, int, int],
        mpp: float | tuple[float, float],
        tile_size: tuple[int, int] = (512, 512),
        pyramid: bool = False,
        colormap: dict[int, str] | None = None,
        compression: TiffCompression | None = TiffCompression.JPEG,
        is_mask: bool = False,
        quality: int | None = 100,
        metadata: dict[str, str] | None = None,
    ):
        """
        Writer based on tifffile.

        Parameters
        ----------
        filename : PathLike
            Filename where to write
        size : tuple
            Size of the image to be written. This is defined as (height, width, num_channels),
            or rather (rows, columns, num_channels) and is important value to get correct.
            In case of a mask with a single channel the value is given by (rows, columns).
        mpp : int, or tuple[int, int]
        tile_size : tuple[int, int]
            Tiff tile_size, defined as (height, width).
        pyramid : bool
            Whether to write a pyramidal image.
        colormap : dict[int, str]
            Colormap to use for the image. This is only used when the image is a mask.
        compression : TiffCompression
            Compressor to use.
        is_mask : bool
            If true a 2x2 maximal filter will be used for the downsampling, otherwise a 2x2 average filter will be used.
        quality : int
            Quality in case a lossy compressor is used.
        metadata : dict[str, str]
            Metadata to write to the tiff file.
        """
        super().__init__(
            filename,
            size,
            mpp,
            tile_size,
            pyramid,
            colormap,
            compression,
            is_mask,
            quality,
            metadata,
        )

    def from_tiles_iterator(self, iterator: Iterator[npt.NDArray[np.int_]]) -> None:
        """
        Generate the tiff from a tiles iterator. The tiles should be in row-major (C-order) order.
        The `dlup.tiling.Grid` class has the possibility to generate such grids using `GridOrder.C`.

        Parameters
        ----------
        iterator : Iterator
            Iterator providing the tiles as numpy arrays.
            They are expected to be (tile_height, tile_width, num_channels) when RGB(A) images or
            (tile_height, tile_width) when masks. The tiles can be smaller at the border.
        """
        filename = pathlib.Path(self._filename)

        native_size = self._size[:-1]
        software = f"dlup {dlup.__version__} (tifffile.py {tifffile.__version__})"
        n_subresolutions = 0
        if self._pyramid:
            n_subresolutions = int(np.ceil(np.log2(np.asarray(native_size) / np.asarray(self._tile_size))).min())
        shapes = [np.floor(np.asarray(native_size) / 2**n).astype(int).tolist() for n in range(0, n_subresolutions + 1)]

        # TODO: add to metadata "axes": "TCYXS", and "SignificantBits": 10,
        metadata = {
            "PhysicalSizeX": self._mpp[0],
            "PhysicalSizeXUnit": "µm",
            "PhysicalSizeY": self._mpp[1],
            "PhysicalSizeYUnit": "µm",
        }
        if self._metadata is not None:
            metadata = {**metadata, **self._metadata}

        # Convert the compression variable to a tifffile supported one.
        _compression = TIFFFILE_COMPRESSION[self._compression.value]

        is_rgb = self._size[-1] in (3, 4)
        if is_rgb and self._colormap is not None:
            raise ValueError("Colormaps only work with integer-valued images (e.g. masks).")
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_filename = pathlib.Path(temp_dir) / filename.name
            tiff_writer = tifffile.TiffWriter(temp_filename, bigtiff=True)
            self._write_page(
                tiff_writer,
                tile_iterator=iterator,
                level=0,
                compression=_compression,
                shapes=shapes,
                is_rgb=is_rgb,
                subifds=None,
                software=software,
                metadata=metadata,
            )

            for level in range(0, n_subresolutions):
                tiff_reader = tifffile.TiffReader(temp_filename)
                page = tiff_reader.pages[level]
                tile_iterator = _tile_iterator_from_page(
                    page,  # type: ignore
                    self._tile_size,
                    shapes[level],
                    is_mask=self._is_mask,
                )
                self._write_page(
                    tiff_writer,
                    tile_iterator=tile_iterator,
                    level=level + 1,
                    compression=_compression,
                    shapes=shapes,
                    is_rgb=is_rgb,
                    subfiletype=1,
                    software=software,
                )
                tiff_reader.close()
            tiff_writer.close()
            shutil.move(str(temp_filename), str(filename))

    def _write_page(
        self,
        tiff_writer: tifffile.TiffWriter,
        tile_iterator: Iterator[npt.NDArray[np.int_]],
        level: int,
        compression: str | None,
        shapes: list[tuple[int, int]],
        is_rgb: bool,
        **options: Any,
    ) -> None:
        native_resolution = 1 / np.array(self._mpp) * 10000
        if self._colormap is not None:
            colorspace = "palette"
        elif is_rgb:
            colorspace = "rgb"
        else:
            colorspace = "minisblack"
        tiff_writer.write(
            tile_iterator,  # noqa
            shape=(*shapes[level], self._size[-1]) if is_rgb else (*shapes[level], 1),
            dtype="uint8",
            resolution=(native_resolution[0] / 2**level, native_resolution[1] / 2**level),
            resolutionunit="CENTIMETER",
            photometric=colorspace,
            compression=compression if not self._quality else (compression, self._quality),  # type: ignore
            tile=self._tile_size,
            colormap=self._colormap,
            **options,
        )


def _tiles_iterator_from_pil_image(
    pil_image: PIL.Image.Image, tile_size: tuple[int, int], order: str | GridOrder = "F"
) -> Generator[npt.NDArray[np.int_], None, None]:
    """
    Given a PIL image return a tile-iterator.

    Parameters
    ----------
    pil_image : PIL.Image
    tile_size : tuple
    order : GridOrder or str

    Yields
    ------
    np.ndarray
        Tile outputted in row-major format
    """
    grid = Grid.from_tiling(
        (0, 0),
        size=pil_image.size[::-1],
        tile_size=tile_size,
        tile_overlap=(0, 0),
        mode=TilingMode.overflow,
        order=order,
    )
    for tile_coordinates in grid:
        arr = np.asarray(pil_image)
        cropped_array = arr[
            tile_coordinates[0] : tile_coordinates[0] + tile_size[0],
            tile_coordinates[1] : tile_coordinates[1] + tile_size[1],
        ].astype("uint8")
        yield cropped_array


def _tile_iterator_from_page(
    page: tifffile.TiffPage,
    tile_size: tuple[int, int],
    region_size: tuple[int, int],
    is_mask: bool = False,
) -> Generator[npt.NDArray[np.int_], None, None]:
    """
    Create an iterator from a tiff page. Useful when writing a pyramidal tiff where the previous page is read to write
    the new page. Each tile will be the downsampled version from the previous version.

    Parameters
    ----------
    page : tifffile.TiffPage
    tile_size : tuple
    region_size : tuple
    is_mask : bool
        Whether the image is a mask (important for the downsampling)

    Yields
    ------
    np.ndarray
        Tile outputted in row-major format
    """
    resized_tile_size = tuple(map(lambda x: x * 2, tile_size))
    grid = Grid.from_tiling(
        (0, 0),
        size=region_size,
        tile_size=resized_tile_size,
        tile_overlap=(0, 0),
        mode=TilingMode.overflow,
        order="F",
    )
    for coordinates in grid:
        # The tile size must be cropped to image bounds
        region_end = coordinates + resized_tile_size
        size = np.clip(region_end, 0, region_size) - coordinates

        vips_tile = get_tile(page, (coordinates[1], coordinates[0]), (size[1], size[0]))

        output = vips_tile.reduce(2, 2, kernel="nearest" if is_mask else "linear").numpy()
        yield output
