"""Recommend qPCR marker genes for any taxon when the user has no specific gene.

Enter an organism or genus (e.g. "Plasmodium", "Haemonchus", "SARS-CoV-2"). This
resolves the lineage via NCBI Taxonomy and maps it -- by domain and phylum, so it works
across the whole tree of life -- to a curated, ranked set of standard qPCR markers, each
annotated with what it's used for, conserved-vs-variable (detect the group vs separate
relatives), copy number (sensitivity) and the conventional region. It also handles
"I want X but NOT relative Y" by setting the SAME marker on the off-target side.

Coverage: viruses (pox + the common families with their standard targets), bacteria,
archaea, apicomplexa, kinetoplastids, Giardia/diplomonads, microsporidia, other protists,
fungi, plants, helminths (nematodes/flatworms/acanths), and other animals -- with a
domain-level fallback (Bacteria->16S, Archaea->16S, Eukaryota->18S, Viruses->polymerase)
so an unusual lineage still lands somewhere sensible. Everything is curated domain
knowledge, not generated. Classification needs the NCBI lookup (network); offline it
falls back to name-keyword matching for common genera, else broad eukaryote markers.
"""
from . import ncbi


def _m(gene, name, q, good_for, copies, level, conservation, chemistry, region=""):
    return dict(gene=gene, name=name, q=q, good_for=good_for, copies=copies,
                level=level, conservation=conservation, chemistry=chemistry, region=region)


# ---------- ribosomal / universal building blocks ----------
def _18S(note=""):
    return _m("18S rRNA", "18S small-subunit ribosomal RNA", "18S ribosomal RNA",
              "Catch the whole group -- the standard sensitive detection target." + (" " + note if note else ""),
              "multicopy \u2192 sensitive", "detection", "conserved", "idt_taqman",
              "short conserved sub-region for qPCR")
def _28S():
    return _m("28S rRNA", "28S large-subunit ribosomal RNA (D1\u2013D3)", "28S ribosomal RNA",
              "Genus / species; ribosomal alternative to 18S/ITS.", "multicopy",
              "identification", "moderate", "idt_taqman", "")
def _ITS(label="ITS", q="internal transcribed spacer", chem="sybr_generic", good="rDNA spacer; the usual species-ID workhorse."):
    return _m(label, "internal transcribed spacer (rDNA)", q, good,
              "multicopy \u2192 sensitive", "identification", "variable", chem, "ITS1 or ITS2 sub-region")
def _COI(chem="sybr_generic"):
    return _m("COI", "cytochrome c oxidase subunit I (mito)", "cytochrome oxidase subunit 1",
              "Mitochondrial barcode; species level (Folmer region).", "mtDNA, multicopy",
              "genotyping", "variable", chem, "Folmer region, ~658 bp")
def _CYTB(chem="sybr_generic"):
    return _m("cytb", "cytochrome b (mitochondrial)", "cytochrome b",
              "Species / lineage discrimination; mitochondrial.", "mtDNA, multicopy",
              "genotyping", "variable", chem, "")

# ---------- apicomplexa ----------
A_18S = _18S(); A_28S = _28S()
A_CYTB = _m("cytb", "cytochrome b (mitochondrial)", "cytochrome b",
            "Separate species / lineages -- the haemosporidian MalAvi barcode.",
            "mtDNA, multicopy", "genotyping", "variable", "parasite_mtdna",
            "avian haemosporidia: ~479 bp Hellgren region (MalAvi)")
A_COI = _COI("parasite_mtdna")
A_COX3 = _m("cox3", "cytochrome c oxidase subunit III (mito)", "cytochrome oxidase subunit 3",
            "Mitochondrial target used in haemosporidian detection assays.",
            "mtDNA, multicopy", "genotyping", "variable", "parasite_mtdna", "")
A_ITS2 = _ITS("ITS2", "internal transcribed spacer 2", "sybr_generic", "Used for Babesia / Theileria typing; variable.")
APICOMPLEXA = [A_18S, A_CYTB, A_COI, A_COX3, A_28S, A_ITS2]

# ---------- kinetoplastids (Trypanosoma, Leishmania) ----------
K_KDNA = _m("kDNA", "kinetoplast minicircle DNA", "kinetoplast minicircle",
            "Kinetoplastid's most sensitive target -- thousands of minicircles per cell.",
            "very high-copy \u2192 very sensitive", "detection", "moderate", "idt_taqman",
            "minicircle conserved region")
K_18S = _18S(); K_GGAPDH = _m("gGAPDH", "glycosomal glyceraldehyde-3-P dehydrogenase", "glycosomal GAPDH",
            "Species / phylogeny in trypanosomatids.", "single-copy", "genotyping", "variable", "idt_taqman", "")
K_ITS = _ITS("ITS", "internal transcribed spacer", "sybr_generic", "Species / typing in Leishmania & trypanosomes.")
KINETO = [K_KDNA, K_18S, K_GGAPDH, K_ITS]

# ---------- Giardia / diplomonads ----------
G_18S = _18S(); G_TPI = _m("tpi", "triosephosphate isomerase", "triosephosphate isomerase",
            "Assemblage genotyping (Giardia).", "single-copy", "genotyping", "variable", "idt_taqman", "")
G_GDH = _m("gdh", "glutamate dehydrogenase", "glutamate dehydrogenase",
           "Assemblage genotyping (Giardia).", "single-copy", "genotyping", "variable", "idt_taqman", "")
G_BG = _m("\u03b2-giardin", "beta-giardin (bg)", "beta-giardin",
          "Assemblage genotyping (Giardia).", "single-copy", "genotyping", "variable", "idt_taqman", "")
GIARDIA = [G_18S, G_TPI, G_GDH, G_BG]

# ---------- microsporidia / other protists ----------
MICRO = [_18S("Microsporidian SSU is the standard detection/ID target."),
         _ITS("ITS", "internal transcribed spacer", "sybr_generic", "Species / strain in microsporidia.")]
PROTIST = [_18S(), _ITS("ITS / ITS2"), _28S(), _COI("sybr_generic")]

# ---------- bacteria ----------
B_16S = _m("16S rRNA", "16S ribosomal RNA", "16S ribosomal RNA",
           "Genus / identification -- the universal bacterial marker.",
           "multicopy \u2192 sensitive", "identification", "conserved", "idt_taqman",
           "V3\u2013V4 or V4 hypervariable region")
B_RPOB = _m("rpoB", "RNA polymerase \u03b2-subunit", "rpoB", "Finer species / strain resolution than 16S.",
            "single-copy", "genotyping", "variable", "idt_taqman", "")
B_GYRB = _m("gyrB", "DNA gyrase subunit B", "gyrB", "Species / strain resolution; single-copy.",
            "single-copy", "genotyping", "variable", "idt_taqman", "")
B_GROEL = _m("groEL (hsp60)", "60-kDa chaperonin", "groEL", "Species ID across many genera.",
             "single-copy", "identification", "moderate", "idt_taqman", "")
B_TUF = _m("tuf", "elongation factor Tu", "tuf", "Species identification; single-copy.",
           "single-copy", "identification", "moderate", "idt_taqman", "")
B_RECA = _m("recA", "recombinase A", "recA", "Phylogeny / species discrimination.",
            "single-copy", "genotyping", "variable", "idt_taqman", "")
B_23S = _m("23S rRNA", "23S ribosomal RNA", "23S ribosomal RNA", "Sensitive detection; ribosomal alternative to 16S.",
           "multicopy \u2192 sensitive", "detection", "conserved", "idt_taqman", "")
B_ITS = _m("16S\u201323S ITS", "16S\u201323S intergenic spacer", "16S 23S intergenic spacer",
           "Strain-level; the most variable rRNA-region marker.", "multicopy",
           "genotyping", "variable", "sybr_generic", "intergenic spacer")
BACTERIA = [B_16S, B_RPOB, B_GYRB, B_GROEL, B_TUF, B_RECA, B_23S, B_ITS]

# ---------- archaea ----------
AR_16S = _m("16S rRNA", "16S ribosomal RNA (archaeal)", "16S ribosomal RNA",
            "Genus / identification -- the universal archaeal marker.", "multicopy \u2192 sensitive",
            "identification", "conserved", "idt_taqman", "archaeal-specific primers")
AR_RPOB = _m("rpoB", "RNA polymerase \u03b2-subunit", "rpoB", "Species / strain resolution.",
             "single-copy", "genotyping", "variable", "idt_taqman", "")
AR_MCRA = _m("mcrA", "methyl-coenzyme M reductase \u03b1", "mcrA",
             "Methanogen-specific functional + phylogenetic marker.", "single-copy",
             "genotyping", "variable", "idt_taqman", "methanogens only")
ARCHAEA = [AR_16S, AR_RPOB, AR_MCRA]

# ---------- fungi ----------
F_ITS = _m("ITS", "internal transcribed spacer (ITS1\u20135.8S\u2013ITS2)", "internal transcribed spacer",
           "The formal fungal barcode; species-level.", "multicopy \u2192 sensitive",
           "identification", "variable", "sybr_generic", "ITS1 or ITS2 sub-region")
F_28S = _m("28S rRNA", "28S large-subunit ribosomal RNA (D1/D2)", "28S ribosomal RNA",
           "Genus / species; complements ITS.", "multicopy", "identification", "moderate", "idt_taqman", "D1/D2 domain")
F_18S = _18S("Broad detection across fungi; conserved.")
F_TEF1 = _m("TEF1\u03b1", "translation elongation factor 1-\u03b1", "translation elongation factor 1 alpha",
            "Species resolution in many genera.", "single-copy", "genotyping", "variable", "idt_taqman", "")
F_BENA = _m("\u03b2-tubulin (benA)", "beta-tubulin", "beta-tubulin",
            "Species ID in Aspergillus / Penicillium.", "single-copy", "genotyping", "variable", "idt_taqman", "")
F_CAM = _m("calmodulin", "calmodulin (CaM)", "calmodulin", "Species ID, esp. Aspergillus / Penicillium.",
           "single-copy", "genotyping", "variable", "idt_taqman", "")
F_RPB2 = _m("RPB2", "RNA polymerase II second-largest subunit", "RPB2",
            "Deep phylogeny / species; single-copy.", "single-copy", "genotyping", "variable", "idt_taqman", "")
FUNGI = [F_ITS, F_28S, F_18S, F_TEF1, F_BENA, F_CAM, F_RPB2]

# ---------- plants ----------
P_RBCL = _m("rbcL", "RuBisCO large subunit (chloroplast)", "rbcL", "Core plant barcode; easy to amplify.",
            "chloroplast, multicopy", "identification", "moderate", "sybr_generic", "")
P_MATK = _m("matK", "maturase K (chloroplast)", "matK", "Plant barcode; more variable than rbcL.",
            "chloroplast, multicopy", "genotyping", "variable", "sybr_generic", "")
P_ITS = _m("ITS / ITS2", "nuclear internal transcribed spacer", "internal transcribed spacer",
           "Nuclear barcode; high resolution.", "multicopy", "genotyping", "variable", "sybr_generic", "")
P_TRNH = _m("trnH-psbA", "trnH-psbA intergenic spacer (chloroplast)", "trnH-psbA",
            "Variable chloroplast spacer for species.", "chloroplast, multicopy", "genotyping",
            "variable", "sybr_generic", "chloroplast spacer")
P_TRNL = _m("trnL (P6)", "trnL intron P6 loop (chloroplast)", "trnL",
            "Short; works on degraded DNA / eDNA / diet.", "chloroplast, multicopy", "detection",
            "moderate", "sybr_generic", "P6 loop, ~10\u2013143 bp")
PLANTS = [P_RBCL, P_MATK, P_ITS, P_TRNH, P_TRNL]

# ---------- helminths (nematodes, flatworms, acanthocephalans) ----------
H_ITS2 = _m("ITS2", "internal transcribed spacer 2 (rDNA)", "internal transcribed spacer 2",
            "The helminth species-ID workhorse; variable rDNA spacer.", "multicopy \u2192 sensitive",
            "identification", "variable", "sybr_generic", "rDNA spacer")
H_ITS1 = _m("ITS1", "internal transcribed spacer 1 (rDNA)", "internal transcribed spacer 1",
            "Species ID; complements ITS2.", "multicopy \u2192 sensitive", "identification", "variable", "sybr_generic", "rDNA spacer")
H_18S = _18S("Broad detection / higher-level phylogeny.")
H_28S = _m("28S rRNA", "28S large-subunit ribosomal RNA (D1\u2013D3)", "28S ribosomal RNA",
           "Genus / family phylogeny.", "multicopy", "identification", "moderate", "idt_taqman", "")
H_COI = _COI("sybr_generic")
H_CYTB = _CYTB("sybr_generic")
H_NAD1 = _m("nad1", "NADH dehydrogenase subunit 1 (mito)", "NADH dehydrogenase subunit 1",
            "Mitochondrial; strong for trematode / cestode species & populations.",
            "mtDNA, multicopy", "genotyping", "variable", "sybr_generic", "")
HELMINTH = [H_ITS2, H_ITS1, H_18S, H_28S, H_COI, H_CYTB, H_NAD1]

# ---------- animals / metazoa (vertebrates, arthropods, molluscs, ...) ----------
Z_COI = _m("COI", "cytochrome c oxidase subunit I (mito)", "cytochrome oxidase subunit 1",
           "The animal DNA barcode; species ID.", "mtDNA, multicopy", "identification", "variable",
           "sybr_generic", "Folmer region, ~658 bp")
Z_CYTB = _CYTB("sybr_generic")
Z_12S = _m("12S rRNA", "12S ribosomal RNA (mitochondrial)", "12S ribosomal RNA",
           "Vertebrate detection / eDNA metabarcoding.", "mtDNA, multicopy", "detection", "moderate",
           "sybr_generic", "MiFish region, ~170 bp")
Z_16S = _m("16S rRNA", "16S ribosomal RNA (mitochondrial)", "16S ribosomal RNA",
           "Species / higher; invertebrate metabarcoding.", "mtDNA, multicopy", "identification",
           "moderate", "sybr_generic", "")
Z_COII = _m("COII", "cytochrome c oxidase subunit II (mito)", "cytochrome oxidase subunit 2",
            "Additional mito marker for species.", "mtDNA, multicopy", "genotyping", "variable", "sybr_generic", "")
Z_DLOOP = _m("D-loop", "mitochondrial control region", "control region",
             "Most variable mtDNA; population-level.", "mtDNA, multicopy", "genotyping", "variable",
             "sybr_generic", "control region")
Z_18S = _18S("Conserved; higher-level detection (esp. invertebrates).")
METAZOA = [Z_COI, Z_CYTB, Z_12S, Z_16S, Z_COII, Z_DLOOP, Z_18S]

# ---------- viruses ----------
V_4B = _m("P4b (fpv167)", "4b core protein gene (Avipoxvirus)", "4b core protein",
          "Standard avipoxvirus detection / clade typing.", "single-copy (large dsDNA virus)",
          "identification", "moderate", "idt_taqman", "fpv167 locus")
POX = [V_4B]


def _vmark(gene, name, q, good, region):
    return _m(gene, name, q, good, "viral genome", "detection", "moderate", "idt_taqman", region)

HERPES = [_vmark("DNA pol (UL30)", "viral DNA polymerase", "DNA polymerase",
                 "Conserved; the usual herpesvirus qPCR target.", "polymerase gene"),
          _vmark("gB", "glycoprotein B", "glycoprotein B", "Type-able structural target.", "glycoprotein B")]
CORONA = [_vmark("RdRp (ORF1ab)", "RNA-dependent RNA polymerase", "RdRp",
                 "Conserved replicase; primary coronavirus target.", "ORF1ab"),
          _vmark("N gene", "nucleocapsid", "nucleocapsid", "Abundant transcript; sensitive screen.", "N gene"),
          _vmark("E gene", "envelope", "envelope", "Conserved screening target (sarbecovirus).", "E gene")]
INFLU = [_vmark("M (matrix)", "matrix gene", "matrix gene",
                "The conserved type-A screening target (CDC assay).", "matrix gene"),
         _vmark("NP", "nucleoprotein", "nucleoprotein", "Subtyping / typing.", "nucleoprotein")]
RETRO = [_vmark("pol", "polymerase (reverse transcriptase)", "pol gene", "Conserved; quantify proviral/RNA load.", "pol"),
         _vmark("gag", "group-specific antigen", "gag gene", "Conserved structural target.", "gag")]
FLAVI = [_vmark("NS5 (RdRp)", "non-structural 5 / RdRp", "NS5", "Conserved replicase; primary flavivirus target.", "NS5"),
         _vmark("5\u02b9UTR", "5\u02b9 untranslated region", "5 UTR", "Highly conserved screening region.", "5\u02b9 UTR")]
ADENO = [_vmark("hexon", "hexon gene", "hexon", "Conserved + type-able adenovirus target.", "hexon")]
PICORNA = [_vmark("5\u02b9UTR", "5\u02b9 untranslated region", "5 UTR", "Conserved enterovirus/rhinovirus screen.", "5\u02b9 UTR"),
           _vmark("3D pol", "3D RNA polymerase", "3D polymerase", "Conserved replicase.", "3D polymerase")]

# taxon-specific gold-standard targets, prepended when the name matches
TOXO_RE = _m("RE (529 bp)", "529-bp repeat element (GenBank AF146527)", "529 bp repeat element",
             "Toxoplasma's most sensitive target.", "high-copy repeat \u2192 very sensitive",
             "detection", "conserved", "idt_taqman", "~529 bp, ~200\u2013300 copies/genome")
TOXO_B1 = _m("B1", "B1 gene", "B1 gene", "Classic multicopy Toxoplasma target.",
             "multicopy (~35 copies) \u2192 sensitive", "detection", "moderate", "idt_taqman", "~35 copies/genome")
MTB_IS6110 = _m("IS6110", "IS6110 insertion sequence", "IS6110",
                "M. tuberculosis complex -- multicopy, very sensitive (complex-specific, not all NTM).",
                "multicopy insertion seq \u2192 very sensitive", "detection", "moderate", "idt_taqman", "MTB-complex specific")
MTB_HSP65 = _m("hsp65", "65-kDa heat-shock protein (groEL2)", "hsp65", "Species ID across mycobacteria.",
               "single-copy", "identification", "moderate", "idt_taqman", "")

_SPECIFIC = [
    (["toxoplasma"], [TOXO_RE, TOXO_B1]),
    (["mycobacterium"], [MTB_IS6110, MTB_HSP65]),
]


def gold_standard_query(organism, gene):
    """If (organism, gene) names a curated taxon-specific gold-standard target -- e.g. Toxoplasma
    'B1' / 'RE', Mycobacterium 'IS6110' / 'hsp65' -- return its marker query so a direct design
    query routes to the free-text marker fetch instead of colliding with an unrelated gene symbol
    (the classic failure: Toxoplasma 'B1' -> human cyclin B1, CCNB1). Returns None when no match,
    so normal protein-coding genes are unaffected."""
    org = (organism or "").lower()
    g = (gene or "").lower().replace("gene", "").replace("-", " ").strip()
    if not g:
        return None
    gtoks = {t for t in g.split() if t}
    for keys, markers in _SPECIFIC:
        if any(k in org for k in keys):
            for m in markers:
                hay = (m["gene"] + " " + m["name"] + " " + m["q"]).lower().replace("-", " ")
                if g in hay or (gtoks & set(hay.split())):
                    return m["q"]
    return None

_RELATIVES = {
    "haemoproteus": ["Plasmodium", "Leucocytozoon"], "plasmodium": ["Haemoproteus", "Leucocytozoon"],
    "leucocytozoon": ["Plasmodium", "Haemoproteus"], "babesia": ["Theileria"], "theileria": ["Babesia"],
}

# ordered most-specific-first; matched against (name + NCBI lineage), lowercased.
_RULES = [
    # --- viruses (pox first, then families) ---
    (["avipox", "fowlpox", "canarypox", "poxvir"], "avipoxvirus / poxvirus", POX,
     "Avipox typing keys on the 4b core protein gene (fpv167)."),
    (["herpesvir", "herpesviridae", "herpesvirales"], "herpesvirus", HERPES,
     "Herpesvirus qPCR usually targets the conserved DNA polymerase; gB for typing."),
    (["coronavir", "coronaviridae", "nidovirales"], "coronavirus", CORONA,
     "RdRp (ORF1ab) is the conserved primary target; N for sensitivity, E for screening."),
    (["orthomyxovir", "influenza"], "influenza virus", INFLU,
     "The matrix (M) gene is the conserved type-A target; HA/NA for subtyping."),
    (["retrovir", "retroviridae", "lentivir", "ortervirales"], "retrovirus", RETRO,
     "pol/gag are conserved; quantify proviral DNA or genomic RNA."),
    (["flavivir", "flaviviridae"], "flavivirus", FLAVI,
     "NS5 (RdRp) and the 5\u02b9UTR are the conserved targets."),
    (["adenovir", "adenoviridae"], "adenovirus", ADENO,
     "The hexon gene is conserved and type-able."),
    (["picornavir", "picornaviridae", "enterovirus", "rhinovirus"], "picornavirus", PICORNA,
     "The 5\u02b9UTR is the conserved enterovirus/rhinovirus screen."),
    # --- cellular: most-specific eukaryote groups before general ones ---
    (["apicomplexa", "haemosporida", "plasmodium", "haemoproteus", "leucocytozoon", "babesia",
      "theileria", "toxoplasma", "eimeria", "cryptosporidium", "cytauxzoon", "sarcocystis",
      "hepatozoon"], "apicomplexan parasite", APICOMPLEXA,
     "18S for sensitive group detection; cytb/COI for species or lineage. For close lineages, "
     "put the diagnostic SNP under the probe -- see the Conservation tab."),
    (["kinetoplast", "trypanosoma", "leishmania", "trypanosomatida", "euglenozoa"], "kinetoplastid parasite",
     KINETO, "kDNA minicircles give the most sensitive detection; gGAPDH/ITS for species/typing."),
    (["giardia", "diplomonad", "hexamita", "spironucleus"], "diplomonad (Giardia-type)", GIARDIA,
     "18S for detection; gdh/tpi/\u03b2-giardin for assemblage genotyping."),
    (["microsporidia", "encephalitozoon", "enterocytozoon", "nosema"], "microsporidian", MICRO,
     "SSU (18S) is the standard detection/ID target; ITS for strain."),
    (["fungi", "ascomycota", "basidiomycota", "saccharomyc", "aspergillus", "candida", "cryptococcus",
      "fusarium", "penicillium", "mucor", "dikarya", "mucoromycota", "chytridiomycota"], "fungus", FUNGI,
     "ITS is the formal fungal barcode; 28S/TEF1\u03b1/benA/CaM resolve where ITS under-resolves."),
    (["viridiplantae", "streptophyta", "magnoliopsida", "embryophyta", "spermatophyta", "chlorophyta"],
     "plant", PLANTS, "rbcL + matK are the core plant barcodes; ITS2/trnH-psbA add resolution; trnL (P6) for degraded DNA."),
    (["platyhelminthes", "trematoda", "cestoda", "monogenea", "nematoda", "acanthocephala",
      "chromadorea", "enoplea", "rhabditida", "strongylida", "ascarid", "spirurida", "trichocephalida",
      "schistosoma", "fasciola", "haemonchus", "ascaris", "trichuris", "echinococcus", "taenia"],
     "helminth (parasitic worm)", HELMINTH,
     "rDNA ITS1/ITS2 is the species-ID workhorse; 18S/28S for phylogeny; mtDNA COI/cytb/nad1 for "
     "species & populations. For a worm vs its host, set the host genus under \u201c\u2026but NOT in\u201d."),
    (["chordata", "vertebrata", "mammalia", "aves", "actinopterygii", "amphibia", "reptilia", "sauropsida",
      "metazoa", "eumetazoa", "arthropoda", "insecta", "mollusca", "cnidaria", "echinodermata", "annelida",
      "porifera"], "animal / metazoan", METAZOA,
     "COI/cytb/12S identify or detect the animal itself. To MEASURE EXPRESSION of a gene, target its "
     "mRNA and normalise to validated reference genes -- see the Ref genes tab."),
    (["alveolata", "stramenopil", "rhizaria", "amoebozoa", "ciliophora", "dinophyceae", "apicomonad",
      "haptista", "cryptophyceae", "discoba", "metamonada", "sar;", "oomycota", "perkinsozoa"],
     "protist", PROTIST, "18S is the universal protist marker for detection/ID; ITS/28S/COI add resolution."),
    (["archaea", "euryarchaeota", "crenarchaeota", "thaumarchaeota", "methanobacteria", "halobacteria",
      "nitrososphaeria", "thermoproteota"], "archaeon", ARCHAEA,
     "16S with archaeal-specific primers; rpoB for finer resolution; mcrA for methanogens."),
    (["; bacteria", "proteobacteria", "firmicutes", "actinobacteria", "bacteroidetes", "spirochaet",
      "cyanobacteria", "mycobacterium", "salmonella", "escherichia", "staphylococc", "streptococc",
      "borrelia", "campylobacter", "listeria", "clostridium", "bacillota", "pseudomonadota"], "bacterium", BACTERIA,
     "16S identifies to genus but often not to species; rpoB/gyrB/groEL resolve finer. For quantification "
     "prefer a single-copy target so copy number is defined."),
]

_GENERIC_EUK = [_18S(), _COI("sybr_generic"), _CYTB("sybr_generic"), _ITS("ITS")]
_GENERIC_NOTE = ("Couldn't place this taxon precisely -- these are broad eukaryote markers. "
                 "A genus or species NCBI recognises (network) gives a sharper recommendation.")
_VIRAL_NOTE = ("This virus has no listed family rule. Target a conserved region of the polymerase "
               "(RdRp / DNA pol) or a family-specific gene from a published assay for it.")

_RANK = {"conserved": 0, "moderate": 1, "variable": 2}

# verbose marker phrases -> the form that actually matches GenBank annotations.
# "Haemoproteus 18S ribosomal RNA" returns 0 records; "Haemoproteus 18S" returns hundreds.
# Applied when building the query string so both the data scan and the Target->assay
# fetch (which the chip fills) find sequences.
_ROBUST = {
    "18S ribosomal RNA": "18S", "28S ribosomal RNA": "28S", "16S ribosomal RNA": "16S",
    "23S ribosomal RNA": "23S", "12S ribosomal RNA": "12S",
    "cytochrome oxidase subunit 1": "cox1", "cytochrome oxidase subunit 2": "cox2",
    "cytochrome oxidase subunit 3": "cox3",
}


def _dedup(markers):
    seen, out = set(), []
    for m in markers:
        if m["gene"] not in seen:
            seen.add(m["gene"]); out.append(m)
    return out


def suggest(name, exclude=None, intent="any"):
    name = (name or "").strip()
    intent = (intent or "any").strip().lower()
    info = ncbi.taxonomy_lineage(name)
    sci, rank, lineage = info if info else (name, "", "")
    hay = (name + " " + lineage).lower()

    group = markers = note = None
    for kws, glabel, mk, gnote in _RULES:
        if any(k in hay for k in kws):
            group, markers, note = glabel, list(mk), gnote
            break

    if markers is None:  # domain-level fallback so every taxon lands somewhere sensible
        if ("virus" in hay) or ("viridae" in hay) or ("viral" in hay):
            group, markers, note = "virus (no universal marker)", [], _VIRAL_NOTE
        elif "archaea" in hay:
            group, markers, note = "archaeon", list(ARCHAEA), _RULES[-1][3]
        elif "bacteria" in hay:
            group, markers, note = "bacterium", list(BACTERIA), "16S for ID; rpoB/gyrB for finer resolution."
        else:
            group, markers, note = "generic eukaryote", list(_GENERIC_EUK), _GENERIC_NOTE

    for kws, extra in _SPECIFIC:
        if any(k in hay for k in kws):
            markers = _dedup(extra + markers)
            break

    if intent in ("detect", "detection", "quantify"):
        markers = sorted(markers, key=lambda m: _RANK.get(m["conservation"], 1))
    elif intent in ("genotype", "genotyping", "discriminate"):
        markers = sorted(markers, key=lambda m: -_RANK.get(m["conservation"], 1))

    base = sci or name
    out = []
    for m in markers:
        qw = _ROBUST.get(m["q"], m["q"])
        out.append(dict(gene=m["gene"], name=m["name"], good_for=m["good_for"], copies=m["copies"],
                        level=m["level"], conservation=m["conservation"], region=m["region"],
                        chemistry=m["chemistry"], qwords=qw, query=("%s %s" % (base, qw)).strip()))

    rels = []
    for key, vals in _RELATIVES.items():
        if key in hay:
            rels = vals; break

    return dict(input=name, resolved=base, rank=rank, lineage=lineage, resolved_ok=bool(info),
                group=group, markers=out, relatives=rels, exclude=(exclude or "").strip(),
                intent=intent, note=note)
