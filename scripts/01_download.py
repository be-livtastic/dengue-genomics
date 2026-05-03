 # Verify Biopython and Pandas installation
try:
    import Bio
    import pandas as pd
    print("Biopython and Pandas are installed correctly.")
except ImportError:
    print("Error: Biopython or Pandas are not installed.")

from __future__ import annotations
py
import csv
import re
import sys
import time
from collections import Counter
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional
from xml.etree import ElementTree as ET

try:
    from Bio import Entrez
except ImportError:
    sys.exit("Biopython not installed. Run: pip install biopython pandas")


# ============================================================================
# CONFIGURATION
# ============================================================================

# Re-use global NCBI_EMAIL and NCBI_API_KEY from previous cells
# NCBI_EMAIL is required by NCBI
# NCBI_API_KEY is optional but recommended
NCBI_EMAIL = "2515050@chester.ac.uk"      # REQUIRED by NCBI
NCBI_API_KEY: Optional[str] = None         # Optional but recommended (value from previous cell)

TAXONOMY_ID = 12637  # Corrected: Dengue virus (verified on NCBI Taxonomy)
INCLUDE_CHROMOSOME = True # include chromosome-level assemblies
SCAFFOLD_FALLBACK_THRESHOLD = 500  # if <500 high-quality, include scaffold-level for more data
MAX_RECORDS = 20000  # Increased safety cap for potentially more Dengue records

# Geographical regions for filtering
TARGET_REGIONS = [
    "caribbean", "central america", "south america",
    "argentina", "belize", "bolivia", "brazil", "chile", "colombia", "costa rica",
    "cuba", "dominica", "dominican republic", "ecuador", "el salvador", "french guiana",
    "grenada", "guadeloupe", "guatemala", "guyana", "haiti", "honduras", "jamaica",
    "martinique", "mexico", "nicaragua", "panama", "paraguay", "peru", "puerto rico",
    "st. lucia", "st. vincent and the grenadines", "suriname", "trinidad and tobago",
    "uruguay", "venezuela"
]

# Output
OUTPUT_DIR_DENGUE = Path("./results")
TIMESTAMP_DENGUE = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H-%M-%SZ")

# ============================================================================
# DATA STRUCTURES
# ============================================================================

@dataclass
class DengueAssemblyRecord:
    assembly_accession: str = ""
    strain_name: str = ""
    assembly_level: str = ""
    biosample: str = ""
    bioproject: str = ""
    isolation_source_original: str = ""
    host: str = ""
    geo_loc_name: str = ""
    country: str = "" # Derived from geo_loc_name
    collection_date: str = ""
    collection_year: str = "" # Derived from collection_date
    submission_date: str = ""
    serotype: str = "" # Extracted from strain_name/notes if possible
    notes: str = ""

# ============================================================================
# NCBI ENTREZ QUERIES (Adapted)
# ============================================================================

def setup_entrez_dengue() -> None:
    if NCBI_EMAIL == "your.email@example.com":
        sys.exit("Set NCBI_EMAIL at the top of this script before running.")
    Entrez.email = NCBI_EMAIL
    Entrez.tool = "DengueFeverScopingReview"
    if NCBI_API_KEY:
        Entrez.api_key = NCBI_API_KEY

def build_search_term_dengue(db: str, include_chromosome: bool, include_scaffold: bool = False, include_geo_filter: bool = True) -> str:
    levels = ['"complete genome"[Assembly Level]']
    if include_chromosome:
        levels.append('"chromosome"[Assembly Level]')
    if include_scaffold:
        levels.append('"scaffold"[Assembly Level]')

    if db == "assembly":
        level_clause = "(" + " OR ".join(levels) + ")"
        base_term = f'txid{TAXONOMY_ID}[Organism:exp] AND {level_clause} AND "latest refseq"[filter] AND all[filter] NOT anomalous[filter]'
    elif db == "nucleotide":
        # For nucleotide, look for 'complete genome' or 'full length' in title or definition
        # Also, filter for RefSeq nucleotide entries specifically (e.g., 'refseq [filter]')
        level_clause = '("complete genome"[Title] OR "full length"[Title] OR "complete sequence"[Title] OR "complete cds"[Title])'
        base_term = f'txid{TAXONOMY_ID}[Organism:exp] AND {level_clause} AND "refseq"[filter]'
    else:
        raise ValueError(f"Unknown database: {db}")

    if include_geo_filter:
        geo_clauses = [f'"{region}"[geo_loc_name]' for region in TARGET_REGIONS]
        geo_term = "(" + " OR ".join(geo_clauses) + ")"
        return f'{base_term} AND {geo_term}'
    else:
        return base_term

def search_ncbi_db(db: str, term: str) -> list[str]:
    print(f"  Searching {db} with Query: {term}")
    handle = Entrez.esearch(db=db, term=term, retmax=MAX_RECORDS)
    result = Entrez.read(handle)
    handle.close()
    time.sleep(0.34 if not NCBI_API_KEY else 0.11) # Add sleep after esearch to respect rate limits
    ids = result["IdList"]
    print(f"  Found {result['Count']} records in {db}, retrieved {len(ids)} UIDs.")
    return ids

def fetch_ncbi_summaries(db: str, uids: list[str], batch_size: int = 200) -> list[dict]:
    all_summaries = []
    for i in range(0, len(uids), batch_size):
        batch = uids[i:i + batch_size]
        print(f"  Fetching summaries from {db} {i + 1}-{i + len(batch)} of {len(uids)}...", end="\r")
        handle = Entrez.esummary(db=db, id=",".join(batch))
        result = Entrez.read(handle)
        handle.close()
        all_summaries.extend(result["DocumentSummarySet"]["DocumentSummary"])
        time.sleep(0.34 if not NCBI_API_KEY else 0.11)  # respect rate limits
    print()
    return all_summaries

def parse_biosample_xml_dengue(xml_str: str) -> dict[str, str]:
    # This function is kept for completeness but its use is deemphasized
    # in main_dengue due to issues with BioSample metadata for viruses.
    fields = {"isolation_source": "", "host": "", "geo_loc_name": "", "collection_date": ""}
    if not xml_str:
        return fields
    try:
        root = ET.fromstring(xml_str)
    except ET.ParseError:
        return fields
    for attr in root.iter("Attribute"):
        name = (attr.get("attribute_name") or attr.get("harmonized_name") or "").lower()
        value = (attr.text or "").strip()
        if not value:
            continue
        if name in {"isolation_source", "isolation source"}:
            fields["isolation_source"] = value
        elif name == "host":
            fields["host"] = value
        elif name in {"geo_loc_name", "geographic location"}:
            fields["geo_loc_name"] = value
        elif name in {"collection_date", "collection date"}:
            fields["collection_date"] = value
    return fields

def fetch_biosample_metadata_dengue(biosample_acc: str) -> dict[str, str]:
    # This function is kept for completeness but its use is deemphasized
    # in main_dengue due to issues with BioSample metadata for viruses.
    if not biosample_acc:
        return {"isolation_source": "", "host": "", "geo_loc_name": "", "collection_date": ""}
    try:
        handle = Entrez.esearch(db="biosample", term=f"{biosample_acc}[accession]")
        result = Entrez.read(handle)
        handle.close()
        if not result["IdList"]:
            return {"isolation_source": "", "host": "", "geo_loc_name": "", "collection_date": ""}
        bs_uid = result["IdList"][0]
        handle = Entrez.efetch(db="biosample", id=bs_uid, rettype="xml")
        xml_str = handle.read().decode("utf-8")
        handle.close()
        return parse_biosample_xml_dengue(xml_str)
    except Exception as e:
        print(f"    [warn] BioSample fetch failed for {biosample_acc}: {e}")
        return {"isolation_source": "", "host": "", "geo_loc_name": "", "collection_date": ""}

def extract_serotype(strain_name: str, additional_text: str = "") -> str:
    # Combine strain_name and any other relevant text for comprehensive search
    text_to_search = f"{strain_name} {additional_text}".lower()

    # Pattern for DENV-1, DENV-2, DENV-3, DENV-4
    serotype_match = re.search(r'denv[ -]?(\d)', text_to_search)
    if serotype_match:
        return f"DENV-{serotype_match.group(1)}"

    # Specific common names (case-insensitive)
    if "dengue virus type 1" in text_to_search: return "DENV-1"
    if "dengue virus type 2" in text_to_search: return "DENV-2"
    if "dengue virus type 3" in text_to_search: return "DENV-3"
    if "dengue virus type 4" in text_to_search: return "DENV-4"
    if "denv1" in text_to_search: return "DENV-1"
    if "denv2" in text_to_search: return "DENV-2"
    if "denv3" in text_to_search: return "DENV-3"
    if "denv4" in text_to_search: return "DENV-4"

    return ""

def summary_to_dengue_record(db: str, summary: dict) -> DengueAssemblyRecord:
    rec = DengueAssemblyRecord()

    if db == "assembly":
        rec.assembly_accession = str(summary.get("AssemblyAccession", "")).strip()
        rec.strain_name = str(summary.get("Organism", "")).strip()

        biosource_info = summary.get("Biosource", {})
        if biosource_info and biosource_info.get("InfraspeciesList"):
            for infra_spec in biosource_info["InfraspeciesList"]:
                if infra_spec.get("Sub_type") == "strain" and infra_spec.get("Sub_value"):
                    if infra_spec["Sub_value"].lower() not in rec.strain_name.lower():
                        rec.strain_name = f"{rec.strain_name} ({infra_spec['Sub_value']})"
                    break

        rec.assembly_level = str(summary.get("AssemblyStatus", "")).strip()
        rec.biosample = str(summary.get("BioSampleAccn", "")).strip()
        rec.bioproject = str(summary.get("BioprojectAccn", "")).strip()
        rec.submission_date = str(summary.get("SubmissionDate", "")).strip()

        if biosource_info:
            isolates = biosource_info.get("Isolates")
            if isolates and isinstance(isolates, list) and len(isolates) > 0:
                first_isolate = isolates[0]
                rec.isolation_source_original = str(first_isolate.get("Source", "")).strip()
                rec.host = str(first_isolate.get("Host", "")).strip()
                rec.geo_loc_name = str(first_isolate.get("GeographicLocation", "")).strip()
                rec.collection_date = str(first_isolate.get("CollectionDate", "")).strip()
            else:
                rec.isolation_source_original = str(biosource_info.get("Source", "")).strip()
                rec.host = str(biosource_info.get("Host", "")).strip()
                rec.geo_loc_name = str(biosource_info.get("GeographicLocation", "")).strip()
                rec.collection_date = str(biosource_info.get("CollectionDate", "")).strip()

    elif db == "nucleotide":
        rec.assembly_accession = str(summary.get("AccessionVersion", "")).strip() # Accession.version for nucleotide
        rec.strain_name = str(summary.get("Title", "")).strip() # Title often contains strain info
        # Nucleotide summaries don't have 'assembly_level' in the same way, can infer from title
        if "complete genome" in rec.strain_name.lower() or "full length" in rec.strain_name.lower():
            rec.assembly_level = "Complete Genome (Nucleotide)"
        else:
            rec.assembly_level = "Nucleotide Sequence"

        rec.biosample = str(summary.get("AcessionVesion","")) # Biosample can be in DBLink
        rec.bioproject = str(summary.get("Bioproject", "")).strip()
        rec.submission_date = str(summary.get("UpdateDate", "")).strip() # Use UpdateDate for nucleotide

        # Attempt to extract common fields from nucleotide summary
        # These fields might be less consistently named than in Assembly Biosource
        rec.isolation_source_original = str(summary.get("Source", "")).strip()
        rec.host = str(summary.get("Host", "")).strip()
        rec.geo_loc_name = str(summary.get("GeographicLocation", "")).strip()
        rec.collection_date = str(summary.get("CreateDate", "")).strip() # Use CreateDate or other date fields

    # Extract collection year
    if rec.collection_date:
        try:
            date_part = rec.collection_date.split('T')[0]
            if re.match(r"^\d{4}-\d{2}-\d{2}$", date_part):
                rec.collection_year = str(datetime.strptime(date_part, '%Y-%m-%d').year)
            elif re.match(r"^\d{4}-\d{2}$", date_part):
                rec.collection_year = str(datetime.strptime(date_part, '%Y-%m').year)
            elif re.match(r"^\d{4}$", date_part):
                rec.collection_year = date_part
        except ValueError:
            pass

    # Extract country from geo_loc_name
    if rec.geo_loc_name:
        country_parts = rec.geo_loc_name.split(':', 1)
        if len(country_parts) > 0:
            rec.country = country_parts[0].strip()

    # Extract serotype
    rec.serotype = extract_serotype(rec.strain_name, rec.isolation_source_original)
    return rec

# ============================================================================
# MAIN PIPELINE (Adapted)
# ============================================================================

def main_dengue() -> None:
    OUTPUT_DIR_DENGUE.mkdir(parents=True, exist_ok=True)
    setup_entrez_dengue()

    print("=" * 70)
    print("Dengue Fever Genome Assembly Curator")
    print(f"Started: {TIMESTAMP_DENGUE}")
    print("=" * 70)

    all_uids = []
    all_summaries = []

    # ----- Step 1.1: Search RefSeq for all global high-quality assemblies (from Assembly DB) -----
    print("\n[1/5] Searching NCBI Assembly database for global high-quality genomes...")
    term_hq_assembly = build_search_term_dengue(db="assembly", include_chromosome=INCLUDE_CHROMOSOME, include_scaffold=False, include_geo_filter=False)
    uids_assembly = search_ncbi_db(db="assembly", term=term_hq_assembly)
    used_term_assembly = term_hq_assembly

    if len(uids_assembly) < SCAFFOLD_FALLBACK_THRESHOLD:
        print(f"\n  Assembly DB high-quality count ({len(uids_assembly)}) below threshold ({SCAFFOLD_FALLBACK_THRESHOLD}). Including global scaffold-level...")
        term_all_assembly = build_search_term_dengue(db="assembly", include_chromosome=True, include_scaffold=True, include_geo_filter=False)
        uids_assembly = search_ncbi_db(db="assembly", term=term_all_assembly)
        used_term_assembly = term_all_assembly

    print(f"  Final Assembly DB selection based on query '{used_term_assembly}' resulted in {len(uids_assembly)} UIDs.")
    all_uids.extend(uids_assembly)

    # ----- Step 1.2: Search RefSeq for all global complete genomes (from Nucleotide DB) -----
    print("\n[2/5] Searching NCBI Nucleotide database for global complete genomes...")
    term_nucleotide = build_search_term_dengue(db="nucleotide", include_chromosome=False, include_scaffold=False, include_geo_filter=False) # Scaffold/chromosome not directly applicable for nucleotide title search
    uids_nucleotide = search_ncbi_db(db="nucleotide", term=term_nucleotide)
    used_term_nucleotide = term_nucleotide
    print(f"  Final Nucleotide DB selection based on query '{used_term_nucleotide}' resulted in {len(uids_nucleotide)} UIDs.")
    all_uids.extend(uids_nucleotide)

    # Remove duplicates from combined UIDs (if any accession is found in both)
    all_uids_set = set(all_uids)
    print(f"\nCombined unique UIDs from Assembly and Nucleotide databases: {len(all_uids_set)}")
    all_uids_list = list(all_uids_set)

    if not all_uids_list:
        print("No assemblies or nucleotide sequences found. Check your query and connection.")
        # Proceed to write empty files if no UIDs were found
        write_dengue_csv([], "dengue_global_genomes.csv")
        write_dengue_metadata_summary([], used_term_assembly + " OR " + used_term_nucleotide, "global_retrieval_metadata.txt", "GLOBAL DENGUE FEVER")
        write_dengue_csv([], "dengue_geo_filtered_genomes.csv")
        write_dengue_metadata_summary([], (used_term_assembly + " OR " + used_term_nucleotide) + " (geo-filtered)", "geo_filtered_retrieval_metadata.txt", "GEO-FILTERED DENGUE FEVER")
        print("\nDone. No records retrieved globally.")
        return # Exit the main function

    # ----- Step 2: Fetch assembly summaries -----
    print(f"\n[3/5] Fetching metadata for {len(uids_assembly)} Assembly records...")
    summaries_assembly = fetch_ncbi_summaries(db="assembly", uids=uids_assembly)

    print(f"\n[4/5] Fetching metadata for {len(uids_nucleotide)} Nucleotide records...")
    summaries_nucleotide = fetch_ncbi_summaries(db="nucleotide", uids=uids_nucleotide)

    records = []
    for s in summaries_assembly:
        rec = summary_to_dengue_record(db="assembly", summary=s)
        if rec.assembly_accession:
            records.append(rec)

    for s in summaries_nucleotide:
        rec = summary_to_dengue_record(db="nucleotide", summary=s)
        if rec.assembly_accession:
            records.append(rec)

    print(f"\n  Finished processing. Kept {len(records)} records for global dataset.")

    # ----- Step 3: no further enrichment from BioSample needed as Biosource is prioritized ----- # Now includes nucleotide metadata parsing
    print(f"\n[5/5] Writing outputs to {OUTPUT_DIR_DENGUE}/")

    # Write global dataset
    if records:
        write_dengue_csv(records, "dengue_global_genomes.csv")
        write_dengue_metadata_summary(records, used_term_assembly + " OR " + used_term_nucleotide, "global_retrieval_metadata.txt", "GLOBAL DENGUE FEVER")
        print(f"  Wrote global dataset: {len(records)} records.")
    else:
        print("  No records available for the global dataset.")
        write_dengue_csv([], "dengue_global_genomes.csv") # Ensure empty file is created
        write_dengue_metadata_summary([], used_term_assembly + " OR " + used_term_nucleotide, "global_retrieval_metadata.txt", "GLOBAL DENGUE FEVER")

    # Apply geographical filter and write separate output
    geo_filtered_records = []
    target_regions_lower = {region.lower() for region in TARGET_REGIONS}

    for rec in records:
        # Check country derived from geo_loc_name
        if rec.country and rec.country.lower() in target_regions_lower:
            geo_filtered_records.append(rec)
        # Also check if any part of geo_loc_name matches a target region
        elif rec.geo_loc_name:
            if any(target in rec.geo_loc_name.lower() for target in target_regions_lower):
                 geo_filtered_records.append(rec)

    if geo_filtered_records:
        write_dengue_csv(geo_filtered_records, "dengue_geo_filtered_genomes.csv")
        write_dengue_metadata_summary(geo_filtered_records, (used_term_assembly + " OR " + used_term_nucleotide) + " (geo-filtered)", "geo_filtered_retrieval_metadata.txt", "GEO-FILTERED DENGUE FEVER")
        print(f"  Wrote geo-filtered dataset: {len(geo_filtered_records)} records.")
    else:
        print("  No records matched the geographical filter criteria for a separate output.")
        write_dengue_csv([], "dengue_geo_filtered_genomes.csv") # Ensure empty file is created
        write_dengue_metadata_summary([], (used_term_assembly + " OR " + used_term_nucleotide) + " (geo-filtered)", "geo_filtered_retrieval_metadata.txt", "GEO-FILTERED DENGUE FEVER")

    print("\nDone.")
    print(f"  Total records processed (global): {len(records)}")
    print(f"  Records matching geographical filter: {len(geo_filtered_records)}")
    print(f"  Output directory:  {OUTPUT_DIR_DENGUE.resolve()}")

def write_dengue_csv(records: list[DengueAssemblyRecord], output_filename: str) -> None:
    out = OUTPUT_DIR_DENGUE / output_filename
    fieldnames = [
        "assembly_accession", "strain_name", "assembly_level",
        "biosample", "bioproject",
        "isolation_source_original", "host", "geo_loc_name", "country",
        "collection_date", "collection_year", "submission_date", "serotype",
        "notes",
    ]
    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for rec in records:
            row = asdict(rec)
            writer.writerow({k: row[k] for k in fieldnames})
    print(f"  Wrote {out.name}")

def write_dengue_metadata_summary(records: list[DengueAssemblyRecord], query: str, output_filename: str, title_prefix: str) -> None:
    out = OUTPUT_DIR_DENGUE / output_filename
    n = len(records) or 1

    levels = Counter(r.assembly_level for r in records)
    biosamples = [r.biosample for r in records if r.biosample]
    bioprojects = [r.bioproject for r in records if r.bioproject]
    duplicate_biosamples = [b for b, c in Counter(biosamples).items() if c > 1]

    # Metadata completeness for prioritized fields
    has_collection_year = sum(1 for r in records if r.collection_year)
    has_country = sum(1 for r in records if r.country)
    has_serotype = sum(1 for r in records if r.serotype)
    has_host = sum(1 for r in records if r.host)

    with out.open("w", encoding="utf-8") as f:
        f.write(f"{title_prefix} RETRIEVAL METADATA & QUALITY CONTROL SUMMARY\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Generated:           {TIMESTAMP_DENGUE}\n")
        f.write(f"Tool:                Entrez/Biopython, taxid {TAXONOMY_ID}\n")
        f.write(f"Query:               {query}\n")
        f.write(f"Records retrieved:   {len(records)}\n\n")

        f.write("ASSEMBLY LEVEL DISTRIBUTION\n")
        f.write("-" * 70 + "\n")
        for lvl, c in levels.most_common():
            f.write(f"  {lvl:<25} {c:>5}  ({c / n * 100:.1f}%)\n")

        f.write("\nPRIORITIZED METADATA COMPLETENESS\n")
        f.write("-" * 70 + "\n")
        f.write(f"  Collection Year      {has_collection_year:>5}/{n}  ({has_collection_year / n * 100:.1f}%)\n")
        f.write(f"  Country              {has_country:>5}/{n}  ({has_country / n * 100:.1f}%)\n")
        f.write(f"  Serotype             {has_serotype:>5}/{n}  ({has_serotype / n * 100:.1f}%)\n")
        f.write(f"  Host                 {has_host:>5}/{n}  ({has_host / n * 100:.1f}%)\n")

        f.write("\nDUPLICATION CHECK\n")
        f.write("-" * 70 + "\n")
        f.write(f"  Unique BioSamples:   {len(set(biosamples))} (of {len(biosamples)})\n")
        f.write(f"  Unique BioProjects:  {len(set(bioprojects))} (of {len(bioprojects)})\n")
        f.write(f"  Duplicate BioSamples: {len(duplicate_biosamples)}\n")
        if duplicate_biosamples:
            f.write("    " + ", ".join(duplicate_biosamples[:20]))
            if len(duplicate_biosamples) > 20:
                f.write(f"  ... and {len(duplicate_biosamples) - 20} more")
            f.write("\n")
    print(f"  Wrote {out.name}")


if __name__ == "__main__":
    main_dengue()
