"""Generate a greyscale .png file of the geometry for each component."""

# python libraries
import json
import os
from PIL import Image, ImageDraw

# 3rd party libraries
import numpy as np
import matplotlib.pyplot as plt

def draw_component_mask(
    scale_to_image_size: bool,
    working_directory: str,
    image_size: int = 512,
    pixels_per_mm: int = 5
) -> None:
    """
    Draw a component mask image from FEMMT geometry data.

    :param scale_to_image_size: If True, scale the output to the specified square
                                image size. If False, preserve the original aspect
                                ratio using the specified pixel density.
    :type scale_to_image_size: bool
    :param working_directory: FEMMT working directory containing the results folder.
    :type working_directory: str
    :param image_size: Width and height of the scaled square image in pixels.
    :type image_size: int
    :param pixels_per_mm: Pixel density used when the original aspect ratio is preserved.
    :type pixels_per_mm: int
    """
    with open(os.path.join(working_directory, "results", "log_coordinates_description.json"), "r") as file:
        data_1 = json.load(file)

    with open(os.path.join(working_directory, "results", "log_electro_magnetic.json"), "r") as file:
        data_2 = json.load(file)

    # Data extraction
    p_outer = [tuple(p) for p in data_1["p_outer"]]
    p_ww = [tuple(p) for p in data_1["p_ww"]]
    thickness_core_insulation = list(data_2["simulation_settings"]["insulation"]["core_insulations"])

    p_air_gap_center = [tuple(p) for p in data_1.get("p_air_gap_center", [])]
    lengths_air_gap = data_1.get("lengths_air_gap", [])

    p_cond_center_1 = [tuple(p[:2]) for p in data_1.get("p_cond_center_1", [])]
    radius_cond_1 = data_1.get("radius_cond_1", 0)

    p_cond_center_2 = [tuple(p[:2]) for p in data_1.get("p_cond_center_2", [])]
    radius_cond_2 = data_1.get("radius_cond_2", 0)

    # Sort points clockwise
    p_outer_sorted = sort_points_clockwise(np.array(p_outer))
    p_ww_sorted = sort_points_clockwise(np.array(p_ww))

    # Determine limits
    all_x = np.concatenate([p_outer_sorted[:, 0], p_ww_sorted[:, 0]])
    all_y = np.concatenate([p_outer_sorted[:, 1], p_ww_sorted[:, 1]])
    min_x, max_x = all_x.min(), all_x.max()
    min_y, max_y = all_y.min(), all_y.max()

    def draw_polygon(coords: np.ndarray, gray_value: int) -> None:
        """
        Draw a filled polygon on the image using grayscale.

        :param coords: List of (x, y) coordinates representing the polygon vertices.
        :type coords: list of tuple
        :param gray_value: Grayscale value to fill the polygon with (0–255).
        :type gray_value: int
        """
        pixels = [to_pixel(x, y) for x, y in coords]
        draw.polygon(pixels, fill=gray_value)

    def to_pixel(x: int, y: int) -> tuple[int, int]:
        """
        Convert real-world coordinates (x, y) to pixel coordinates for image rendering.

        :param x: X-coordinate in real-world units.
        :type x: int
        :param y: Y-coordinate in real-world units.
        :type y: int
        :return: Corresponding pixel coordinates as (x_px, y_px).
        :rtype: tuple[int, int]
        """
        x_px = int((x - min_x) * scale_x)
        y_px = int((y - min_y) * scale_y)
        return x_px, height_px - y_px

    if scale_to_image_size:
        width_px = image_size
        height_px = image_size
        scale_x = image_size / (max_x - min_x)
        scale_y = image_size / (max_y - min_y)

        scales = {
            "scale_x": scale_x,
            "scale_y": scale_y
        }

        with open(os.path.join(working_directory, "scales.json"), "w") as f:
            json.dump(scales, f, indent=4)

        save_path = os.path.join(working_directory, "component_mask", "scaled", "component_mask.png")

        if not os.path.exists(save_path):
            os.makedirs(save_path)

    else:
        pixels_per_unit = pixels_per_mm * 1000  # mm → m
        width_px = int((max_x - min_x) * pixels_per_unit)
        height_px = int((max_y - min_y) * pixels_per_unit)
        scale_x = scale_y = pixels_per_unit

        save_path = os.path.join(working_directory, "component_mask", "unscaled", "component_mask.png")

        if not os.path.exists(save_path):
            os.makedirs(save_path)

    mask = Image.new("L", (width_px, height_px), 0)
    draw = ImageDraw.Draw(mask)

    # 1. Core
    draw_polygon(p_outer_sorted, 50)

    # 2. Core Insulation
    (x0, y0), (x1, y1), (x2, y2), (x3, y3) = p_ww
    t_top, t_bottom, t_left, t_right = thickness_core_insulation
    core_insulation = np.array([
        (x2 + t_left, y2 - t_top), (x3 - t_right, y2 - t_top),
        (x3 - t_right, y1 + t_bottom), (x2 + t_left, y1 + t_bottom)
    ])

    draw_polygon(core_insulation, 100)

    # 3. Air Gaps
    for center, length in zip(p_air_gap_center, lengths_air_gap, strict=False):
        x, y = center
        half_length = length / 2

        right_bound = p_ww[0][0]
        left_bound = 0

        px_right, _ = to_pixel(right_bound, 0)
        px_left, _ = to_pixel(left_bound, 0)
        _, py_top = to_pixel(0, y - half_length)
        _, py_bottom = to_pixel(0, y + half_length)

        draw.polygon([
            (px_right, py_top), (px_right, py_bottom),
            (px_left, py_bottom), (px_left, py_top)
        ], fill=150)

    # 4. Winding window
    draw_polygon(p_ww_sorted, 0)

    # 5. Conductors type 1
    for cx, cy in p_cond_center_1:
        cx_px, cy_px = to_pixel(cx, cy)
        rx = int(radius_cond_1 * scale_x)
        ry = int(radius_cond_1 * scale_y)
        draw.ellipse([cx_px - rx, cy_px - ry, cx_px + rx, cy_px + ry], fill=200)

    # 6. Conductors type 2
    for cx, cy in p_cond_center_2:
        cx_px, cy_px = to_pixel(cx, cy)
        rx = int(radius_cond_2 * scale_x)
        ry = int(radius_cond_2 * scale_y)
        draw.ellipse([cx_px - rx, cy_px - ry, cx_px + rx, cy_px + ry], fill=255)

    # Show and save
    plt.imshow(mask, cmap='gray')
    plt.axis('off')
    # plt.show()

    mask.save(save_path)

    img_array = np.array(mask, dtype=np.uint8)
    csv_path = save_path.replace(".png", ".csv")
    np.savetxt(csv_path, img_array, fmt="%d", delimiter=";")

def sort_points_clockwise(pts: np.ndarray) -> np.ndarray:
    """
    Sort a set of 2D points in clockwise order around their centroid.

    :param pts: Array of 2D points with shape (N, 2).
    :type pts: np.ndarray
    :return: Points sorted in clockwise order.
    :rtype: np.ndarray
    """
    center = np.mean(pts, axis=0)
    angles = np.arctan2(pts[:, 1] - center[1], pts[:, 0] - center[0])
    return pts[np.argsort(angles)]
