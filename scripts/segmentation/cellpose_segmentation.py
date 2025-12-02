#
# Cell segmentation using cellpose on morphology images from 10X Xenium data
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage:
#        python cellpose_segmentation.py

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0 Import libraries
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import tifffile
import shapely.geometry as geometry
import zarr
from cellpose import models, io, plot
from tifffile import imread
from shapely.geometry import Polygon
from shapely.affinity import translate, scale
from shapely.errors import TopologicalError
from rasterio import features
from sklearn.metrics import jaccard_score
from skimage.measure import regionprops_table

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Load the data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set up plotting aesthetics and variables
sns.set(style='whitegrid')
slide = 'slide_6'
os.chdir(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT_new/data/raw/{slide}/')

# Loading and inspecting the morphology image
image_data_path = 'morphology.ome.tif'
image = tifffile.imread(image_data_path)
print(f"Image shape: {image.shape}")

# Inspect additional metadata regarding the image
with tifffile.TiffFile(image_data_path) as tif:
    ome_metadata = tif.ome_metadata
    print(ome_metadata)

# Preview the transcripts data
data = pd.read_parquet('transcripts.parquet')
print(data.head())

# Select dimensions to analyze
print(f"Image dimensions: {image.ndim}")
if image.ndim == 5:
    # Example shape: (Time, Z, Channels, Height, Width)
    # Select the first time point, z-slice, and channel
    image_channel = image[0, 0, 0, :, :]
elif image.ndim == 4:
    # Example shape: (Z, Channels, Height, Width)
    image_channel = image[0, 0, :, :]
elif image.ndim == 3:
    # Example shape: (Channels, Height, Width)
    image_channel = image[0, :, :]
else:
    # Already a 2D image
    image_channel = image

print(f'Image shape: {image.shape}')
print(f'Image channel shape: {image_channel.shape}')
print(f'np.max: {np.max(image, axis=0)}')

image_channel = np.max(image, axis=0)
image_channel = image_channel.astype(np.uint16)

plt.figure(figsize=(8, 8))
plt.imshow(image_channel)
plt.title('Selected Image for Segmentation')
plt.axis('off')
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT_new/segmentation/plots/{slide}_image.png', dpi=300)

# Use the whole image for segmentation (no ROI)
# image_channel already contains the 2D image we'll use
image_for_seg = image_channel
image_height, image_width = image_channel.shape

# Use scaling factor from the image metadata (pixels per µm)
x_scale = 1 / 0.2125  # pixels per µm (approx 4.7059 px/µm)
y_scale = x_scale

# Map transcripts to the whole image (convert µm -> pixels)
# Keep original `data` but add pixel coordinates for mapping
data['x_pixel'] = (data['x_location'] * x_scale).astype(int)
data['y_pixel'] = (data['y_location'] * y_scale).astype(int)

print(f"Using full image for segmentation: shape={image_for_seg.shape}")

plt.figure(figsize=(8, 8))
plt.imshow(image_for_seg)
plt.title('Whole Image for Segmentation')
plt.axis('off')
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT_new/segmentation/plots/{slide}_whole_image_for_segmentation.png', dpi=300)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Cell pose
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Select pretrained model
model = models.CellposeModel(model_type='cyto3')    # add gpu=False to run on CPU
cell_diameter_um = 10  # diameter of a cell in µm
cell_diameter_pixels = cell_diameter_um * x_scale
print(f"Estimated cell diameter in pixels: {cell_diameter_pixels}")

# Run segmentation
# - `masks`: a labeled mask array where each detected cell has a unique identifier,
# - `flows`: which provides information about cell boundary flows,  
# - `styles`: representing style vectors for detected objects, and
# - `diams`: the diameter used in the model (helpful if it has been automatically adjusted).
masks, flows, styles = model.eval(image_for_seg, diameter=cell_diameter_pixels,)
print(f"Number of cells detected in image: {masks.max()}")

# Locate transcripts to the corresponding cells (use pixel coordinates added to `data`)
x_indices = data['x_pixel'].values
y_indices = data['y_pixel'].values

# Image dimensions
img_height, img_width = image_for_seg.shape

# Clip indices to image bounds
x_indices = np.clip(x_indices, 0, img_width - 1)
y_indices = np.clip(y_indices, 0, img_height - 1)

# Assign each transcript to a segmented cell based on its pixel coordinates
cell_labels = masks[y_indices, x_indices]
data['cellpose_cell_id'] = cell_labels
print(data.head())

# Keep only transcripts assigned to a cell
assigned_data = data[data['cellpose_cell_id'] > 0].copy()
print(f"Total transcripts in image: {len(data)}")
print(f"Assigned transcripts: {len(assigned_data)}")

# Visualize the segmentation overlayed on the image
fig = plt.figure(figsize=(50, 50))
plot.show_segmentation(fig, image_for_seg, masks, flows[0])
plt.title('Cellpose Segmentation', fontsize=40)
fig.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT_new/segmentation/plots/{slide}_segmentation.png', dpi=300)

# Group by cell and gene to get expression counts
expression_per_cell = assigned_data.groupby(['cellpose_cell_id', 'feature_name']).size().reset_index(name='count')
expression_matrix = expression_per_cell.pivot(index='cellpose_cell_id', columns='feature_name', values='count').fillna(0)
print(expression_matrix.head())

# Visualize the transcript locations overlayed on the image
plt.figure(figsize=(20, 20))
plt.imshow(image_for_seg, cmap='gray')
plt.scatter(
    x_indices[data['cellpose_cell_id'] > 0],
    y_indices[data['cellpose_cell_id'] > 0],
    c='red', s=5, label='Transcripts'
)
plt.title('Transcripts Mapped to Segmented Cells')
plt.axis('off')
plt.legend()
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT_new/segmentation/plots/{slide}_transcript_loc.png', dpi=300)

# Visualize GOI
gene_of_interest = 'CD3D'
if gene_of_interest in expression_matrix.columns:
    cell_ids = expression_matrix.index.values
    expression_values = expression_matrix[gene_of_interest].values

    #get centroids of cells
    from skimage.measure import regionprops
    properties = regionprops(masks)
    centroids = np.array([prop.centroid for prop in properties])
    cell_labels = np.array([prop.label for prop in properties])

    #create a mapping from cell label to centroid
    centroid_dict = {label: centroid for label, centroid in zip(cell_labels, centroids)}

    #get centroids for the cells in expression_matrix
    cell_centroids = np.array([centroid_dict.get(cell_id, (np.nan, np.nan)) for cell_id in cell_ids])

    plt.figure(figsize=(8, 8))
    plt.imshow(roi_image, cmap='gray')
    plt.scatter(
        cell_centroids[:, 1],  # x-coordinates
        cell_centroids[:, 0],  # y-coordinates
        c=expression_values,
        cmap='viridis',
        s=5,
        edgecolors='k',
        label=f'Expression of {gene_of_interest}'
    )
    plt.title(f'Expression of {gene_of_interest}')
    plt.axis('off')
    plt.colorbar(label='Expression Level')
    plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT_new/segmentation/plots/{slide}_test_ROI_{gene_of_interest}.png', dpi=300)
else:
    print(f"{gene_of_interest} not found in expression matrix.")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Compare cellpose with 10X segmentation using Jaccard index
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cells_data_path = 'cells.csv.gz'
cells_data = pd.read_csv(cells_data_path)
print(cells_data.head())

boundaries_path = 'nucleus_boundaries.csv.gz'  
boundaries = pd.read_csv(boundaries_path)
print(boundaries.head())

# Use full image to compare 10X segmentation (no ROI)
# Convert image bounds to micrometers
image_height, image_width = image_for_seg.shape
x_start_um = 0.0
y_start_um = 0.0
x_end_um = image_width / x_scale
y_end_um = image_height / y_scale

# Filter cells whose centroids are within the ROI 
cells_in_image = cells_data[
    (cells_data['x_centroid'] >= x_start_um) &
    (cells_data['x_centroid'] < x_end_um) &
    (cells_data['y_centroid'] >= y_start_um) &
    (cells_data['y_centroid'] < y_end_um)].copy()
print(f"Number of cells in ROI from original segmentation: {len(cells_in_image)}")

# Filter cell boundaries for cells in ROI
boundaries_in_image = boundaries[boundaries['cell_id'].isin(cells_in_image['cell_id'])].copy()
cell_polygons = {}
for cell_id, group in boundaries_in_image.groupby('cell_id'):
    x_coords = group['vertex_x'].values
    y_coords = group['vertex_y'].values
    coords = list(zip(x_coords, y_coords))
    try:
        polygon = Polygon(coords)
        if not polygon.is_valid:
            # Attempt to fix invalid polygons
            polygon = polygon.buffer(0)
        cell_polygons[cell_id] = polygon
    except TopologicalError as e:
        print(f"Could not create polygon for cell {cell_id}: {e}")

# Function to transform geometries to pixel coordinates
def geometry_to_pixel_coords(geometry):
    # scale geometry from micrometers to pixels
    geometry_scaled = scale(geometry, xfact=x_scale, yfact=y_scale, origin=(0, 0))
    return geometry_scaled

# Apply transformation to all cell polygons
cell_polygons_px = {cell_id: geometry_to_pixel_coords(geom) for cell_id, geom in cell_polygons.items()}

# Map cell_id strings to integer labels
cell_id_to_label = {cell_id: idx+1 for idx, cell_id in enumerate(cell_polygons_px.keys())}
label_to_cell_id = {idx+1: cell_id for idx, cell_id in enumerate(cell_polygons_px.keys())}

# Prepare shapes for rasterization
shapes = [(geom, cell_id_to_label[cell_id]) for cell_id, geom in cell_polygons_px.items()]

# Create an empty mask and rasterize the shapes
original_masks = np.zeros_like(image_for_seg, dtype=np.uint16)
original_masks = features.rasterize(
    shapes,
    out_shape=original_masks.shape,
    fill=0,
    all_touched=True,
    dtype=np.uint16)

# Lets compare the cellpose segmentation with 10x segmentation side by side
plt.figure(figsize=(16, 8))

#Cellpose 
plt.subplot(1, 2, 1)
plt.imshow(roi_image, cmap='gray')
plt.imshow(masks, alpha=0.5, cmap='jet')
plt.title('Cellpose Segmentation')
plt.axis('off')

#Original 
plt.subplot(1, 2, 2)
plt.imshow(roi_image, cmap='gray')
plt.imshow(original_masks, alpha=0.5, cmap='jet')
plt.title('Original 10x Segmentation')
plt.axis('off')

plt.tight_layout()
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT_new/segmentation/plots/{slide}_test_seg_10X_vs_cellpose.png',dpi=300)


# Overlay the cellpose segmentation on the 10x segmentation to see how well they align
plt.figure(figsize=(8, 8))
plt.imshow(roi_image, cmap='gray')
plt.imshow((original_masks > 0).astype(int), cmap='Blues', alpha=0.5, label='Original')
plt.imshow((masks > 0).astype(int), cmap='Reds', alpha=0.5, label='Cellpose')
plt.title('Overlay of Segmentations')
plt.axis('off')
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT_new/segmentation/plots/{slide}_test_seg_10X_vs_cellpose_overlay.png', dpi=300)

# Convert masks to binary masks (cells vs background)
cellpose_mask_binary = (masks > 0).astype(int)
original_mask_binary = (original_masks > 0).astype(int)

# Flatten the masks for metric computation
cellpose_mask_flat = cellpose_mask_binary.flatten()
original_mask_flat = original_mask_binary.flatten()

jaccard = jaccard_score(original_mask_flat, cellpose_mask_flat)
print(f'Jaccard Index: {jaccard:.4f}')

def dice_coefficient(y_true, y_pred):
    intersection = np.sum(y_true * y_pred)
    sum_union = np.sum(y_true) + np.sum(y_pred)
    dice = 2 * intersection / sum_union
    return dice

dice = dice_coefficient(original_mask_flat, cellpose_mask_flat)
print(f'Dice Coefficient: {dice:.4f}')

cellpose_cell_count = masks.max()
original_cell_count = original_masks.max()

print(f"Number of cells detected by Cellpose: {cellpose_cell_count}")
print(f"Number of cells in original segmentation: {original_cell_count}")

#cellpose cell areas
cellpose_props = regionprops_table(masks, properties=['area'])
cellpose_areas = pd.DataFrame(cellpose_props)
cellpose_areas['method'] = 'Cellpose'

#original cell areas
original_props = regionprops_table(original_masks, properties=['area'])
original_areas = pd.DataFrame(original_props)
original_areas['method'] = 'Original'

#combine data
areas_df = pd.concat([cellpose_areas, original_areas], ignore_index=True)

# compare cell sizes between the two segmentation methods
plt.figure(figsize=(10, 6))
sns.kdeplot(data=areas_df, x='area', hue='method', common_norm=False)
plt.xlabel('Cell Area (pixels)')
plt.ylabel('Density')
plt.title('Cell Size Distribution Comparison')
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT_new/segmentation/plots/{slide}test_seg_10X_vs_cellpose_cellarea.png',dpi=300)


# Save the results to a new file
# tifffile.imwrite(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT_new/segmentation/segmentation_results/{slide}_cellpose_masks_roi.tif', masks.astype(np.uint16))
# tifffile.imwrite(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT_new/segmentation/segmentation_results/{slide}_original_masks_roi.tif', original_masks.astype(np.uint16))

# metrics = pd.DataFrame({
#     'Metric': ['Jaccard Index', 'Dice Coefficient'],
#     'Value': [jaccard, dice]
# })
# metrics.to_csv(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT_new/segmentation/segmentation_results/{slide}_segmentation_comparison_metrics.csv', index=False)

# areas_df.to_csv(f'/data/spatial_workshop/day1/Xenium_V1_FFPE_TgCRND8_17_9_months_outs/{slide}_cell_size_comparison.csv', index=False)
