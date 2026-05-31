"""Recommend qPCR marker genes for a taxon when the user has no specific gene.

Given just an organism or genus (e.g. "Plasmodium"), resolve its lineage via NCBI
Taxonomy and map it to a curated, ranked set of standard qPCR marker genes — the
kind a methods reviewer would expect — each with what it is good for, its copy
number (which drives sensitivity), a suggested chemistry profile, and a ready-to-run
Target->assay query string. Falls back to keyword matching on the raw name when the
taxonomy lookup is unavailable (offline), so common genera still work.

The recommendations are curated domain knowledge, not guesses: rRNA cistron markers
(18S/16S/ITS/28S) are multicopy and sensitive for detection/identification;
mitochondrial markers (cytb/COI/12S) resolve species and are multicopy; single-copy
housekeeping genes (rpoB/gyrB/TEF1) give finer strain/species resolution.
"""
from . import ncbi

# ---- marker templates (q = the words to put in the NCBI query) ----
M_18S  = dict(gene="18S rRNA", name="18S small-subunit ribosomal RNA", q="18S ribosomal RNA",
              good_for="Genus-level detection; the standard sensitive target.",
              copies="multicopy \u2192 sensitive", level="genus / detection", chemistry="idt_taqman")
M_CYTB = dict(gene="cytb", name="cytochrome b (mitochondrial)", q="cytochrome b",
              good_for="Species / lineage discrimination (e.g. MalAvi haemosporidian lineages).",
              copies="mtDNA, multicopy", level="species / genotyping", chemistry="parasite_mtdna")
M_COI  = dict(gene="COI", name="cytochrome c oxidase subunit I (mitochondrial)", q="cytochrome oxidase subunit 1",
              good_for="Species-level barcode (Folmer region).",
              copies="mtDNA, multicopy", level="species / barcode", chemistry="parasite_mtdna")
M_16S  = dict(gene="16S rRNA", name="16S ribosomal RNA", q="16S ribosomal RNA",
              good_for="Genus / identification; the universal bacterial marker.",
              copies="multicopy \u2192 sensitive", level="genus / detection", chemistry="idt_taqman")
M_RPOB = dict(gene="rpoB", name="RNA polymerase \u03b2-subunit", q="rpoB",
              good_for="Finer species / strain resolution than 16S.",
              copies="single-copy", level="species", chemistry="idt_taqman")
M_GYRB = dict(gene="gyrB", name="DNA gyrase subunit B", q="gyrB",
              good_for="Species / strain resolution; single-copy.",
              copies="single-copy", level="species", chemistry="idt_taqman")
M_ITS  = dict(gene="ITS", name="internal transcribed spacer (ITS1\u20135.8S\u2013ITS2)", q="internal transcribed spacer",
              good_for="The formal fungal barcode; species-level.",
              copies="multicopy \u2192 sensitive", level="species / barcode", chemistry="sybr_generic")
M_28S  = dict(gene="28S rRNA", name="28S large-subunit ribosomal RNA (D1/D2)", q="28S ribosomal RNA",
              good_for="Genus / species; complements ITS.",
              copies="multicopy", level="genus / species", chemistry="idt_taqman")
M_TEF1 = dict(gene="TEF1\u03b1", name="translation elongation factor 1-\u03b1", q="translation elongation factor 1 alpha",
              good_for="Species resolution in many fungal genera.",
              copies="single-copy", level="species", chemistry="idt_taqman")
M_12S  = dict(gene="12S rRNA", name="12S ribosomal RNA (mitochondrial)", q="12S ribosomal RNA",
              good_for="Vertebrate detection / eDNA metabarcoding (e.g. MiFish).",
              copies="mtDNA, multicopy", level="genus / detection", chemistry="sybr_generic")
M_RBCL = dict(gene="rbcL", name="RuBisCO large subunit (chloroplast)", q="rbcL",
              good_for="Core plant barcode.",
              copies="chloroplast, multicopy", level="genus / species", chemistry="sybr_generic")
M_MATK = dict(gene="matK", name="maturase K (chloroplast)", q="matK",
              good_for="Plant barcode; more variable than rbcL.",
              copies="chloroplast, multicopy", level="species", chemistry="sybr_generic")
M_4B   = dict(gene="P4b (fpv167)", name="4b core protein gene (Avipoxvirus)", q="4b core protein",
              good_for="Standard avipoxvirus detection / clade typing.",
              copies="single-copy (large dsDNA virus)", level="genus / clade", chemistry="idt_taqman")

# ordered most-specific-first; matched against (name + NCBI lineage), lowercased
_RULES = [
    (["avipox", "fowlpox", "canarypox", "poxvir"], "avipoxvirus / poxvirus", [M_4B],
     "Avipox typing keys on the 4b core protein gene (the fpv167 locus). Other viruses "
     "have no universal marker \u2014 target a conserved region of the polymerase or a "
     "family-specific gene from a published assay."),
    (["apicomplexa", "haemosporida", "plasmodium", "haemoproteus", "leucocytozoon", "babesia",
      "theileria", "toxoplasma", "eimeria", "cryptosporidium", "cytauxzoon", "sarcocystis",
      "hepatozoon"], "apicomplexan parasite", [M_18S, M_CYTB, M_COI],
     "18S for sensitive genus-level detection; cytb/COI for species or lineage. For tightly "
     "related lineages, put the diagnostic SNP under the probe \u2014 use the Conservation tab "
     "to find positions conserved in your target but divergent from off-targets."),
    (["fungi", "ascomycota", "basidiomycota", "saccharomyc", "aspergillus", "candida",
      "cryptococcus", "fusarium", "penicillium", "mucor", "dikarya"], "fungus", [M_ITS, M_28S, M_TEF1],
     "ITS is the formal fungal barcode; 28S (D1/D2) and TEF1\u03b1 help where ITS under-resolves."),
    (["bacteria", "proteobacteria", "firmicutes", "actinobacteria", "bacteroidetes", "spirochaet",
      "mycobacterium", "salmonella", "escherichia", "staphylococc", "streptococc", "borrelia",
      "campylobacter", "listeria", "clostridium"], "bacterium", [M_16S, M_RPOB, M_GYRB],
     "16S identifies to genus but often not to species; rpoB/gyrB resolve finer. For "
     "quantification prefer a single-copy target so copy number is defined."),
    (["viridiplantae", "streptophyta", "magnoliopsida", "embryophyta", "spermatophyta"], "plant",
     [M_RBCL, M_MATK, M_ITS],
     "rbcL + matK are the core plant barcodes; ITS2 adds resolution."),
    (["chordata", "vertebrata", "mammalia", "aves", "actinopterygii", "amphibia", "reptilia",
      "metazoa", "arthropoda", "mollusca", "nematoda", "insecta"], "animal / metazoan",
     [M_COI, M_CYTB, M_12S],
     "COI/cytb/12S are for identifying or detecting the animal itself. If instead you want to "
     "MEASURE EXPRESSION of a gene, target that gene's mRNA and normalise to validated "
     "reference genes \u2014 see the Ref genes tab."),
]
_GENERIC = [M_18S, M_COI, M_CYTB, M_ITS]
_GENERIC_NOTE = ("Couldn't place this taxon precisely, so these are broad eukaryote markers. "
                 "Give a genus or a name NCBI recognises for a sharper recommendation.")


def _viral_other(hay):
    return ("virus" in hay or "viridae" in hay or "viral" in hay)


def suggest(name):
    name = (name or "").strip()
    info = ncbi.taxonomy_lineage(name)
    sci, rank, lineage = info if info else (name, "", "")
    hay = (name + " " + lineage).lower()

    group, markers, note = None, None, None
    for kws, glabel, mk, gnote in _RULES:
        if any(k in hay for k in kws):
            group, markers, note = glabel, mk, gnote
            break
    if markers is None and _viral_other(hay):
        group = "virus (no universal marker)"
        markers = []
        note = ("Viruses have no universal qPCR marker. Target a conserved region of the "
                "polymerase (e.g. RdRp / DNA pol) or a family-specific gene taken from a "
                "published assay for your virus.")
    if markers is None:
        group, markers, note = "generic eukaryote", _GENERIC, _GENERIC_NOTE

    base = sci or name
    out = []
    for m in markers:
        out.append(dict(gene=m["gene"], name=m["name"], good_for=m["good_for"],
                        copies=m["copies"], level=m["level"], chemistry=m["chemistry"],
                        query=("%s %s" % (base, m["q"])).strip()))
    return dict(input=name, resolved=base, rank=rank, lineage=lineage,
                resolved_ok=bool(info), group=group, markers=out, note=note)
