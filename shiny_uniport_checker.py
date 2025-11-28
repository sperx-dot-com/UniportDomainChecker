from shiny import App, ui, render, reactive
import requests

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

#UI Frontend
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

#Server Backend
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

app = App(app_ui, server)

if __name__ == "__main__":
    from shiny import run_app
    run_app(app)
