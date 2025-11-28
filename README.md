# UniProt Feature Lookup Tool

A simple Shiny for Python web application that queries the UniProt API to retrieve protein feature annotations overlapping a specified amino acid sequence region.

## Features

- Query UniProt entries by accession ID
- Specify sequence regions (start/end positions)
- View overlapping protein features (domains, binding sites, PTMs, etc.)
- Interactive table view with detailed feature information
- Clean, responsive interface

## Installation

```bash
pip install shiny requests
```

## Usage

Run the application:

```bash
shiny run uniprot_app.py
```

Then open your browser to `http://127.0.0.1:8000`

## Example

Try the default values:
- **UniProt ID:** Q7Z4F1
- **Start position:** 494
- **End position:** 500

Click "Lookup Features" to see overlapping annotations.

## How It Works

The app fetches protein data from the [UniProt REST API](https://www.uniprot.org/help/api) and filters features that overlap with your specified sequence region. Chain features (representing the full protein) are excluded from results.

## Requirements

- Python 3.7+
- shiny
- requests

## License

MIT

## Data Source

All protein data is retrieved from [UniProt](https://www.uniprot.org/), licensed under CC BY 4.0.