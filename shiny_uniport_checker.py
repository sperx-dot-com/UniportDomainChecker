from shiny import App, ui, render, reactive
import requests
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from io import BytesIO
import base64

# =============================================================================
# Core Functionality: UniProt Feature Retrieval
# =============================================================================

def get_feature_descriptions(uniprot_id: str, start: int, end: int) -> str:
    """
    Return comma separated UniProt feature descriptions overlapping a region.
    
    Args:
        uniprot_id: UniProt accession ID (e.g., 'Q7Z4F1')
        start: Start position of the region (1-based, inclusive)
        end: End position of the region (1-based, inclusive)
    
    Returns:
        Comma-separated string of feature descriptions with types in parentheses,
        or 'Unknown' if no features found
    
    Raises:
        ValueError: If the UniProt entry cannot be fetched
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    r = requests.get(url, timeout=10)
    if not r.ok:
        raise ValueError(f"Failed to fetch UniProt entry: {uniprot_id} (Status: {r.status_code})")
    data = r.json()
    descriptions = []
    for feature in data.get("features", []):
        ftype = feature.get("type", "").strip()
        if ftype.lower() == "chain":
            continue  # skip global protein name
        loc = feature.get("location", {})
        fstart = int(loc.get("start", {}).get("value", -1))
        fend = int(loc.get("end", {}).get("value", -1))
        # Check overlap with user region
        if not (fend < start or fstart > end):
            desc = feature.get("description", "").strip()
            label = desc or ftype
            formatted = f"{label} ({ftype})"
            if formatted not in descriptions:
                descriptions.append(formatted)
    return ", ".join(descriptions) if descriptions else "Unknown"


def get_features_detailed(uniprot_id: str, start: int, end: int) -> dict:
    """
    Return detailed structured data about UniProt features overlapping a region.
    
    Args:
        uniprot_id: UniProt accession ID (e.g., 'Q7Z4F1')
        start: Start position of the region (1-based, inclusive)
        end: End position of the region (1-based, inclusive)
    
    Returns:
        Dictionary containing:
            - uniprot_id: The queried UniProt ID
            - protein_name: Full protein name
            - region: String representation of the queried region
            - features: List of dicts with keys: type, description, start, end
    
    Raises:
        ValueError: If the UniProt entry cannot be fetched
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    r = requests.get(url, timeout=10)
    if not r.ok:
        raise ValueError(f"Failed to fetch UniProt entry: {uniprot_id} (Status: {r.status_code})")
    
    data = r.json()
    overlapping_features = []
    
    for feature in data.get("features", []):
        ftype = feature.get("type", "").strip()
        if ftype.lower() == "chain":
            continue
        
        loc = feature.get("location", {})
        fstart = int(loc.get("start", {}).get("value", -1))
        fend = int(loc.get("end", {}).get("value", -1))
        
        # Check overlap with user region
        if not (fend < start or fstart > end):
            overlapping_features.append({
                "type": ftype,
                "description": feature.get("description", "").strip(),
                "start": fstart,
                "end": fend
            })
    
    # Extract protein name
    protein_name = "Unknown"
    try:
        protein_name = data.get("proteinDescription", {}).get(
            "recommendedName", {}
        ).get("fullName", {}).get("value", "Unknown")
    except (KeyError, AttributeError):
        pass
    
    return {
        "uniprot_id": uniprot_id,
        "protein_name": protein_name,
        "region": f"{start}-{end}",
        "features": overlapping_features
    }

# =============================================================================
# Shiny App
# =============================================================================

app_ui = ui.page_fluid(
    ui.panel_title("UniProt Feature Lookup Tool"),
    ui.layout_sidebar(
        ui.sidebar(
            ui.input_text(
                "uniprot_id",
                "UniProt ID",
                value="Q7Z4F1",
                placeholder="e.g., Q8L999, P12345"
            ),
            ui.input_numeric(
                "start_pos",
                "Start position (1-based)",
                value=494,
                min=1
            ),
            ui.input_numeric(
                "end_pos",
                "End position (inclusive)",
                value=500,
                min=1
            ),
            ui.input_action_button(
                "run", 
                "Lookup Features",
                class_="btn-primary w-100"
            ),
            ui.hr(),
            ui.p(
                "Enter a UniProt accession ID and a sequence region to find overlapping features.",
                class_="text-muted small"
            ),
            width=350
        ),
        ui.navset_card_tab(
            ui.nav_panel(
                "Results",
                ui.output_text_verbatim("result", placeholder=True)
            ),
            ui.nav_panel(
                "Feature Details",
                ui.output_ui("details")
            ),
            ui.nav_panel(
                "Visualization",
                ui.output_ui("plot")
            ),
            ui.nav_panel(
                "About",
                ui.markdown("""
                ### About This Tool
                
                This application queries the UniProt API to retrieve protein feature annotations 
                that overlap with a specified amino acid sequence region.
                
                **Features included:**
                - Domain annotations
                - Binding sites
                - Post-translational modifications
                - Secondary structure elements
                - And more...
                
                **Note:** Chain features (representing the full protein) are excluded from results.
                
                **Data source:** [UniProt](https://www.uniprot.org/)
                """)
            )
        )
    )
)

def server(input, output, session):
    # Reactive value to store detailed results
    feature_data = reactive.Value(None)
    
    @reactive.effect
    @reactive.event(input.run)
    def fetch_features():
        try:
            start = int(input.start_pos())
            end = int(input.end_pos())
            uniprot_id = input.uniprot_id().strip()
            
            if not uniprot_id:
                feature_data.set({"error": "Please enter a UniProt ID."})
                return
            
            if end < start:
                feature_data.set({"error": "End position must be ≥ start position."})
                return
            
            # Use the detailed function to get structured data
            data = get_features_detailed(uniprot_id, start, end)
            feature_data.set(data)
            
        except Exception as e:
            feature_data.set({"error": str(e)})
    
    @output
    @render.text
    def result():
        data = feature_data.get()
        if data is None:
            return "Click 'Lookup Features' to begin..."
        
        if "error" in data:
            return f"Error: {data['error']}"
        
        if not data["features"]:
            return f"No features found overlapping region {data['region']} in {data['uniprot_id']}"
        
        # Generate simple comma-separated output from detailed data
        descriptions = []
        for f in data["features"]:
            label = f["description"] or f["type"]
            formatted = f"{label} ({f['type']})"
            if formatted not in descriptions:
                descriptions.append(formatted)
        
        return ", ".join(descriptions)
    
    @output
    @render.ui
    def details():
        data = feature_data.get()
        if data is None:
            return ui.p("No data loaded yet.", class_="text-muted")
        
        if "error" in data:
            return ui.div(
                ui.p(f"Error: {data['error']}", class_="text-danger")
            )
        
        if not data["features"]:
            return ui.div(
                ui.h5(f"Protein: {data['protein_name']}"),
                ui.p(f"No features found in region {data['region']}", class_="text-muted")
            )
        
        # Create a formatted table of features
        rows = []
        for f in data["features"]:
            rows.append(
                ui.tags.tr(
                    ui.tags.td(f["type"], style="font-weight: 500;"),
                    ui.tags.td(f["description"] or "—"),
                    ui.tags.td(f"{f['start']}-{f['end']}", style="text-align: right;")
                )
            )
        
        return ui.div(
            ui.h5(f"Protein: {data['protein_name']}"),
            ui.p(f"Region queried: {data['region']}", class_="text-muted"),
            ui.hr(),
            ui.h6(f"Found {len(data['features'])} overlapping feature(s):"),
            ui.tags.table(
                ui.tags.thead(
                    ui.tags.tr(
                        ui.tags.th("Type"),
                        ui.tags.th("Description"),
                        ui.tags.th("Position", style="text-align: right;")
                    )
                ),
                ui.tags.tbody(*rows),
                class_="table table-sm table-striped"
            )
        )

    @output
    @render.ui
    def plot():
        data = feature_data.get()
        if data is None:
            return ui.p("No data loaded yet.", class_="text-muted")
        
        if "error" in data:
            return ui.div(
                ui.p(f"Error: {data['error']}", class_="text-danger")
            )
        
        if not data["features"]:
            return ui.div(
                ui.p(f"No features to visualize in region {data['region']}", class_="text-muted")
            )
        
        # Create the plot
        fig, ax = plt.subplots(figsize=(12, max(6, len(data["features"]) * 0.5)))
        
        # Define colors for different feature types
        color_map = {
            "Domain": "#3498db",
            "Repeat": "#9b59b6",
            "Region": "#e74c3c",
            "Site": "#f39c12",
            "Binding site": "#e67e22",
            "Modified residue": "#1abc9c",
            "Helix": "#16a085",
            "Turn": "#27ae60",
            "Beta strand": "#2ecc71",
            "Topological domain": "#34495e",
            "Transmembrane": "#95a5a6",
            "Signal peptide": "#d35400",
            "Propeptide": "#c0392b",
        }
        default_color = "#7f8c8d"
        
        # Get protein sequence length (estimate from max feature end)
        max_pos = max([f["end"] for f in data["features"]] + [int(data["region"].split("-")[1])])
        
        # Plot each feature as a horizontal bar
        y_pos = 0
        labels = []
        colors_used = {}
        
        for feature in sorted(data["features"], key=lambda x: x["start"]):
            ftype = feature["type"]
            color = color_map.get(ftype, default_color)
            colors_used[ftype] = color
            
            # Draw feature bar
            start = feature["start"]
            end = feature["end"]
            width = end - start + 1
            
            ax.barh(y_pos, width, left=start, height=0.8, 
                   color=color, edgecolor='black', linewidth=0.5, alpha=0.7)
            
            # Label
            label = feature["description"] if feature["description"] else ftype
            if len(label) > 30:
                label = label[:27] + "..."
            labels.append(f"{label} ({start}-{end})")
            
            y_pos += 1
        
        # Highlight queried region
        region_parts = data["region"].split("-")
        query_start = int(region_parts[0])
        query_end = int(region_parts[1])
        ax.axvspan(query_start, query_end, alpha=0.15, color='red', zorder=-1)
        
        # Formatting
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels)
        ax.set_xlabel("Amino Acid Position", fontsize=12, fontweight='bold')
        ax.set_ylabel("Features", fontsize=12, fontweight='bold')
        ax.set_title(f"{data['protein_name']} - Feature Map", fontsize=14, fontweight='bold', pad=20)
        ax.set_xlim(0, max_pos + 10)
        ax.grid(axis='x', alpha=0.3, linestyle='--')
        
        # Add legend for feature types
        legend_patches = [mpatches.Patch(color=color, label=ftype, alpha=0.7) 
                         for ftype, color in sorted(colors_used.items())]
        legend_patches.append(mpatches.Patch(color='red', alpha=0.15, label='Queried Region'))
        ax.legend(handles=legend_patches, loc='center left', bbox_to_anchor=(1, 0.5), 
                 frameon=True, fancybox=True, shadow=True)
        
        plt.tight_layout()
        
        # Convert plot to image
        buf = BytesIO()
        plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode()
        plt.close(fig)
        
        return ui.div(
            ui.h5(f"Feature Map: {data['protein_name']}"),
            ui.p(f"Queried region: {data['region']} (highlighted in red)", class_="text-muted"),
            ui.HTML(f'<img src="data:image/png;base64,{img_base64}" style="max-width: 100%; height: auto;">'),
            style="padding: 10px;"
        )

app = App(app_ui, server)

if __name__ == "__main__":
    from shiny import run_app
    run_app(app)
