# Ion Channel Visualizer

A comprehensive tool for analyzing and visualizing ion channel models from ModelDB and ICG (Ion Channel Genealogy) databases.

## Overview

This project provides tools to:
- Download and process ion channel models from ModelDB
- Process and analyze ICG (Ion Channel Genealogy) data
- Calculate similarity metrics between models using Levenshtein ratios
- Create network representations of ion channel relationships
- Generate visualizations and analysis figures
- Provide an interactive web-based visualizer

## Project Structure

- `download_modelDB_models/` - Scripts for downloading and processing ModelDB models
- `process_ICG_data/` - ICG data processing and analysis
- `process_modelDB_ion_channel_models/` - ModelDB specific processing
- `compute_levenshtein_ratio/` - Similarity calculation tools
- `create_network_representation/` - Network graph generation
- `combine_modelDB_wth_ICG_data/` - Data integration utilities
- `thesis_figures/` - Analysis and visualization scripts
- `visualizer/` - Original D3.js visualization implementation
- `nextjs-app/` - Modern Next.js web application with interactive visualizer
- `supermodels/` - Supermodel analysis and markdown generation
- `ICG/` - ICG channel data files

## Features

- **Data Collection**: Automated downloading and processing of ion channel models from ModelDB
- **Similarity Analysis**: Levenshtein ratio calculations for model comparison
- **Network Analysis**: Graph-based representation of ion channel relationships
- **Visualization**: Interactive web interface and static figure generation
- **Data Integration**: Unified comprehensive dictionary combining ModelDB and ICG data

## Interactive Visualizer

The project includes a modern web-based visualizer built with Next.js that provides interactive exploration of ion channel networks.

### Getting Started

1. **Navigate to the Next.js app:**
   ```bash
   cd nextjs-app
   ```

2. **Install dependencies:**
   ```bash
   npm install
   ```

3. **Start the development server:**
   ```bash
   npm run dev
   ```

4. **Open the visualizer:**
   Visit `http://localhost:5000/visualizer` in your browser

### Visualization Features

#### **Network Display**
- **Force-directed graph** showing ion channel models as nodes and relationships as links
- **Dynamic node sizing** based on number of identical models
- **Group splitting** that organizes nodes by filter combinations when multiple filters are active
- **Responsive design** with mobile detection (shows desktop-only message on small screens)

#### **Interactive Filtering**
- **Ion Channel Classes**: Filter by K, Na, Ca, IH, KCa, or Other channels
- **Similarity Score**: Adjust threshold (0-100%) to show/hide links based on code similarity
- **Supermodel Filters**: Toggle Supermodel 1 and Supermodel 2 entries
- **ICG Entry Toggle**: Show/hide ICG (Ion Channel Genealogy) entries

#### **Node Interaction & Selection**

- **Regular Click**: Select/deselect individual nodes
  - Selected nodes appear with blue highlight and dark border
  - Click again to deselect

- **Ctrl/Cmd+Click**: Toggle source nodes
  - Source nodes appear in **gold color** with dark border
  - Use for manual selection of comparison sources

- **Shift+Click**: Smart subgraph selection
  - Automatically selects entire connected component (subgraph)
  - Identifies the **chronologically earliest node** as the source (gold)
  - All other connected nodes become targets (blue with border)
  - Uses publication year and ModelDB ID to determine the original model

#### **Source Code Comparison**

The visualizer includes powerful source code analysis features:

##### **Keyboard Shortcuts**
- **`D` key**: Generate side-by-side code diff after selecting nodes
- **`←` / `→` arrow keys**: Navigate between multiple diff comparisons
- **`Esc` key**: Close the diff viewer

##### **Code Fetching**
- Automatically downloads `.mod` files from ModelDB S3 storage
- Caches downloaded files to avoid repeated requests
- Uses original `unique_modelDB_mod_id` from node metadata

##### **Diff Generation**
- Creates HTML side-by-side comparisons highlighting:
  - **Green background**: Added lines
  - **Red background**: Deleted lines  
  - **Yellow background**: Modified lines
  - **White background**: Unchanged lines
- Professional modal interface with navigation controls
- Shows current diff index (e.g., "1 / 3" for multiple comparisons)

#### **Advanced Workflows**

1. **Compare Individual Models**:
   - Ctrl+Click multiple nodes to mark as sources (gold)
   - Click other nodes to select as targets (blue)
   - Press `D` to generate diffs

2. **Analyze Model Evolution**:
   - Shift+Click any node to auto-select its evolutionary tree
   - The earliest model becomes source, derivatives become targets
   - Press `D` to see how the model evolved over time

3. **Filter-Based Analysis**:
   - Use filters to isolate specific ion channel types
   - Multiple active filters create visual node groups
   - Each group gets positioned using optimal circle packing

#### **Technical Features**
- **Responsive Layout**: Full viewport usage on desktop/tablet, mobile warning
- **Local Diff API**: Server-side diff generation (no external dependencies)
- **Dynamic SVG Sizing**: Visualization adapts to container dimensions
- **Window Resize Handling**: Automatically adjusts when browser is resized

### Data Requirements

Place the following JSON files in the `public/` directory:
- `network_data.json` - Node and link data
- `fixed_location_circles.json` - Optimal group positioning data

### Browser Compatibility

- Modern browsers supporting ES6+ and SVG
- Desktop/tablet recommended (768px+ width)
- Mobile users see responsive warning message

## License

[License information to be added]

## Contact

[Contact information to be added]