#!/usr/bin/env python3
r"""How the input library was normalised before running this script, and the
shape of the -k alias file it expects.

#### Download library. ####
# cat lib.fa | sed 's/DNA transposon/DNA_transposon/g' | sed 's/EnSpm\/CACTA/EnSpm_CACTA/g' | sed 's/LTR Retrotransposon/LTR_retrotransposon/g' > lib2.fa
# cat lib2.fa | sed 's/Non-LTR Retrotransposon/non_LTR_retrotransposon/g' | sed 's/Transposable Element/repeat_fragment/g' > lib3.fa 
# cat lib3.fa | sed 's/Non-LTR_retrotransposon/non_LTR_retrotransposon/g' | sed 's/LTR retrotransposon/Class_I/g' | sed 's/Endogenous Retrovirus/Endogenous_Retrovirus_LTR_retrotransposon/g' | sed 's/endogenous Retrovirus/Endogenous_Retrovirus_LTR_retrotransposon/g' | sed 's/Penelope\/Poseidon/Penelope_retrotransposon/g' |  sed 's/Nonautonomous/repeat_fragment/g' |  sed 's/Nonautonomous/repeat_fragment/g' > lib4.fa
# cat lib4.fa | sed 's/ERV-2_PM-I/Endogenous_Retrovirus_LTR_retrotransposon/g' |  sed 's/endogenous retrovirus/Endogenous_Retrovirus_LTR_retrotransposon/g' | sed 's/Inverted repeat/terminal_inverted_repeat/g' | sed 's/Repetitive element/repeat_region/g' > lib5.fa
# cat lib5.fa | sed 's/endogenous transposon/Endogenous_Retrovirus_LTR_retrotransposon/g' |  sed 's/Gypsy-like endogenous retorvirus/Endogenous_Retrovirus_LTR_retrotransposon/g' | sed 's/Telomeric region repeat/telomeric_repeat/g' | sed 's/Repeat region/repeat_region/g' | awk -F'\t' 'BEGIN{OFS="\t"} /^>/ && $1 ~ /ERV/ { $2="Endogenous_Retrovirus_LTR_retrotransposon" } { print }'  > lib6.fa
# cat lib6.fa | awk -F'\t' 'BEGIN{OFS="\t"} /^>/ && $1 ~ /L1/ { $2="L1_LINE_retrotransposon" } { print }' > lib7.fa

# Run as 'python fasta_to_RepeatMasker.py -i lib7.fa -o lib_repeatmasker.fa -k SO.txt -l log --unmatched-to-unknown'

# SO.txt looks like this:
#######################
####### Contents ######
#######################
## Sequence_Ontology	SO_ID	Alias
centromeric_repeat	SO:0001797	centromeric_repeat,Cent,Cent/CentC,CentC,CentM,Centro/tandem,Cent/centromeric_repeat
knob	SO:0002257	knob,knob/knob180,knob/TR-1
satellite_DNA	SO:0000005	satellite_DNA,Satellite/rice,satellite,minisatellite,microsatellite,Satellite/Satellite,Satellite/Y-chromosome,Satellite,Satellite/Y_chromosome,SAT,MSAT
telomeric_repeat	SO:0001496	telomeric_repeat,telomere,telomeric,telomere/telomere
subtelomere	SO:0001997	subtelomere,subtelomere/4-12-1
low_complexity	SO:0001004	low_complexity,Low_complexity,low_complexity_region,Simple_repeat,Simple_repeat/NA
chloroplast_DNA	SO:0001033	chloroplast_DNA,chloroplast/chloroplast
mitochondrial_DNA	SO:0001032	mitochondrial_DNA,mitochondrion/mitochondrion

## higher order
repeat_region	SO:0000657	repeat_region
repeat_fragment	SO:0001050	repeat_fragment,Unknown,unknown,unknown/unknown,NA/NA,Unknown/NA,Unknown/unknown,Unspecified,repeat/unknown,repeat/Unknown,Unknown/Unknown,unknown/NA
retrotransposon	SO:0000180	Class_I,RNA_transposon,retrotransposon,Retroposon,Retroposon/L1_dep
DNA_transposon	SO:0000182	Class_II,DNA_transposon,DNA/unknown,DNA/Unknown,DNA,DNA/NA
snRNA	SO:0000274	snRNA,snRNA/NA
scRNA	SO:0000013	scRNA

## rDNA and NOR
rRNA_gene	SO:0002360	rRNA_gene,rDNA/45S,rRNA,rDNA
rDNA_intergenic_spacer_element	SO:0001860	rDNA_intergenic_spacer_element,rDNA/spacer,rDNA/IGS
cytosolic_2S_rRNA	SO:0002337	cytosolic_2S_rRNA,rDNA/2s_rRNA,rDNA/2S,2S_rRNA,2s_rRNA,rDNA/2s
cytosolic_5S_rRNA	SO:0000652	cytosolic_5S_rRNA,rDNA/5s_rRNA,rDNA/5S,5S_rRNA,5s_rRNA,rDNA/5s
cytosolic_5_8S_rRNA     SO:0000375      cytosolic_5_8S_rRNA,rDNA/5.8s_rRNA,rDNA/5.8S_rRNA,5.8s_rRNA,5_8s_rRNA,5_8S_rRNA_gene,rDNA/5.8S,rDNA/5.8s,5.8S_rRNA,rDNA/5_8S,5_8S_rRNA
cytosolic_16S_rRNA	SO:0001000	cytosolic_16S_rRNA,rDNA/16s_rRNA,rDNA/16S_rRNA,rDNA/16S,rDNA/16s,16S_rRNA,16s_rRNA
cytosolic_18S_rRNA      SO:0000407      cytosolic_18S_rRNA,rDNA/18s_rRNA,rDNA/18S_rRNA,18s_rRNA,rDNA/18S,rDNA/18s,18S_rRNA
cytosolic_23S_rRNA	SO:0001001	cytosolic_23S_rRNA,rDNA/23s_rRNA,rDNA/23S_rRNA,rDNA/23S,rDNA/23s,23S_rRNA,23s_rRNA
cytosolic_25S_rRNA      SO:0001002      cytosolic_25S_rRNA,rDNA/25s_rRNA,rDNA/25S_rRNA,25S_rRNA,rDNA/25S,rDNA/25s,25s_rRNA
cytosolic_28S_rRNA	SO:0000653	cytosolic_28S_rRNA,rDNA/28s_rRNA,rDNA/28S_rRNA,rDNA/28S,rDNA/28s,28S_rRNA,28s_rRNA

## TIR DNA transposons
terminal_inverted_repeat_element	SO:0000208	terminal_inverted_repeat_element,TIR/unknown,DNA/DTX,DTX,DNAauto/unknown,DNAnona/unknown,TIR/Unknown,DNAauto/Unknown,DNAnona/Unknown
MITE	SO:0000338	MITE,TIR/MITE,MITE/unknown,MITE/Unknown
CACTA_TIR_transposon	SO:0002285	CACTA_TIR_transposon,DNAauto/CACTA,DNAnona/CACTA,DNAauto/CACTG,DNAnona/CACTG,DNA/DTC,MITE/DTC,DTC,dSpm,CACTC,En-Spm,EnSpm,CMC-EnSpm,DNA/CACTA,DNA/CACTG,MITE/CACTA,MITE/CACTG,TIR/EnSpm_CACTA,DNA/EnSpm_CACTA,DNA/CMC-EnSpm,DNA/CMC_Chapaev_3,DNA/CMC_EnSpm,DNA/CMC-Chapaever,DNA/CMC-Chapaev-3,DNA/CMC-EnSpmo,DNA/CMC-Miragen
hAT_TIR_transposon	SO:0002279	hAT_TIR_transposon,DNAauto/hAT,DNAnona/hAT,MITE/DTA,DNA/DTA,DTA,hAT,Ac-Ds,Ac/Ds,hAT-Ac,DNA/hAT,MITE/hAT,TIR/hAT,DNA/hAT-Ac,DNA/hAT-Tag1,DNA/hAT-Blackjack,DNA/hAT-Charlie,DNA/hAT-hAT5,DNA/hAT-hAT6,DNA/hAT-hobo,DNA/hAT-Tip100,DNA/hAT_Ac,DNA/hAT_Blackjack,DNA/hAT_Charlie,DNA/hAT_hAT5,DNA/hAT_hAT6,DNA/hAT_hobo,DNA/hAT_Pegasus,DNA/hAT_Tag1,DNA/hAT_Tip100,DNA/hAT-hATm,DNA/hAT-AcTigger,DNA/hAT_AcTigger,DNA/hAT-hAT1,DNA/hAT_hAT1,DNA/hAT-hATm-hAT-hybrid,DNA/hAT-hATm-hAT_hybrid,DNA/hAT-hATw,DNA/hAT_hATw,DNA/hAT-hATx,DNA/hAT_hATx,DNA/hAT-Pegasus,DNA/hAT-Restless,DNA/hAT_Restless,DNA/hAT-Tag11,DNA/hAT_Tag11,DNA/hAT-Tip100?,DNA/hAT_Tip100?,DNA/hAT-Tol2,DNA/hAT_Tol2
Mutator_TIR_transposon	SO:0002280	Mutator_TIR_transposon,DNAauto/MULE,DNAnona/MULE,DNAnona/MULEtir,MITE/DTM,DNA/DTM,DTM,Mutator,MuDR,DNA/MULE,DNA/MULEtir,MITE/MULE,MITE/MULEtir,MULEtir,TIR/MuDR_Mutator,DNA/Mutator,DNA/MuDR,DNA/MULE-MuDR,DNA/MULE-NOF,DNA/MULE_MuDR,DNA/MULE_NOF,DNA/MuLE-NOF?,DNA/MuLE_NOF?,DNA/MuLE-MuDR,DNA/MuLE_MuDR,DNA/MuLE-F,DNA/MuLE_F,DNA/MUDR,DNA/MULE-MuDR?
PIF_Harbinger_TIR_transposon	SO:0002284	PIF_Harbinger_TIR_transposon,DNAnona/Tourist,MITE/Tourist,MITE/DTH,DNA/DTH,DTH,PIF-Harbinger,PIF/Harbinger,Harbinger,Tourist,DNA/Tourist,TIR/PIF_Harbinger,DNA/Harbinger,DNA/PIF-Harbinger,DNA/PIF,DNA/PIF-ISL2EU,DNA/PIF_ISL2EU,DNA/PIF_Harbinger
Tc1_Mariner_TIR_transposon	SO:0002278	Tc1_Mariner_TIR_transposon,stowaway,Stowaway,DNA/DTT,MITE/Stow,MITE/DTT,DTT,Tc1-Mariner,Tc1_Mariner,Tc1_mariner,Tc1/Mariner,TcMar-Stowaway,DNAauto/MLE,DNAnona/MLE,DNA/MLE,MITE/MLE,TIR/Tc1_Mariner,DNA/Tc1-Mariner,DNA/TcMar-Pogo,DNA/Mariner,DNA/TcMar-Stowaway,DNA/TcMar-Mariner,DNA/TcMar,DNA/TcMar-Fot1,DNA/TcMar-ISRm11,DNA/TcMar-Tc1,DNA/TcMar-Tigger,DNA/TcMar_Fot1,DNA/TcMar_ISRm11,DNA/TcMar_Mariner,DNA/TcMar_Tc1,DNA/TcMar_Tc2,DNA/TcMar_Tigger,DNA/TcMar_Pogo,DNA/Tc1_Mariner,DNA/TcMar_Stowaway,DNA/TcMareSL,DNA/TcMar-Ant1,DNA/TcMar_Ant1,DNA/TcMar-Gizmo,DNA/TcMar_Gizmo,DNA/TcMar-IS630,DNA/TcMar_IS630,DNA/TcMar-IS885,DNA/TcMar_IS885,DNA/TcMar-m44,DNA/TcMar_m44,DNA/TcMar-Mogwai,DNA/TcMar_Mogwai,DNA/TcMar-Sagan,DNA/TcMar_Sagan,DNA/TcMar-Tc2,DNA/TcMar-Tc4,DNA/TcMar_Tc4,Mariner/Tc1,DNA/Mariner-Tc1,DNA/MarinerTc1
P_TIR_transposon	SO:0001535	P_TIR_transposon,P-element,P_element,DNA/DTP,TIR/P,DNA/P
piggyBac_TIR_transposon	SO:0002283	piggyBac_TIR_transposon,PiggyBac,DNA/DTB,MITE/DTB,TIR/PiggyBac,TIR/piggyBac,piggyBac,DNA/PiggyBac,DNA/PiggyBacen,DNA/piggyBac
polinton	SO:0001170	polinton,maverick,Maverick,DNA/Maverick,TIR/Maverick,Polinton
Transib_TIR_transposon	SO:0002282	Transib_TIR_transposon,transib,DNA/DTR,MITE/DTR,DNA/CMC-Transib,DNA/CMC-Transibway,DNA/Transib
Merlin_TIR_transposon	SO:0002281	Merlin_TIR_transposon,Merlin,DNA/DTE,MITE/DTE,TIR/Merlin,DNA/Merlin,DNA/Merlin1
terminal_inverted_repeat	SO:0000481	terminal_inverted_repeat,TIR

## nonTIR DNA transposons
Crypton_YR_transposon	SO:0002277	Crypton_YR_transposon,Crypton,DNA/DYC,DYC,DNA/Crypton,DNA/CryptonA,DNA/Crypton_A,DNA/Crypton-A,DNA/CryptonI,DNA/Crypton_I,DNA/Crypton-I,DNA/CryptonS,DNA/Crypton_S,DNA/Crypton-S,DNA/CryptonV,DNA/Crypton_V,DNA/Crypton-V,DNA/CryptonF,DNA/Crypton_F,DNA/Crypton-F,DNA/CryptonH,DNA/Crypton_H,DNA/Crypton-H
helitron	SO:0000544	helitron,DNAauto/Helitron,DNAnona/Helitron,DNA/Helitron,Helitron,RC/Helitron,Unknown/Helitron-2,nonDNA/Helitron,Helitron/NA

## LTR retrotransposons
LTR_retrotransposon	SO:0000186	LTR_retrotransposon,LTR/unknown,LTR/Solo,LTR/Unknown,LTR/solo,LTR,LTR/NA
Retrovirus_LTR_retrotransposon	SO:0002267	Retrovirus_LTR_retrotransposon,LTR/retrovirus,retrovirus,LTR/RLR,RLR,LTR/Retrovirus
TRIM	SO:0002261	TRIM,LTR/TRIM
LARD	SO:0002260	LARD,LTR/LARD
Copia_LTR_retrotransposon	SO:0002264	Copia_LTR_retrotransposon,LTR/Copia,LTR/RLC,RLC,Copia,Ty1,LTR/Ty1,LTR/Copia?
Gypsy_LTR_retrotransposon	SO:0002265	Gypsy_LTR_retrotransposon,LTR/Gypsy,LTR/RLG,RLG,Gypsy,Ty3,LTR/CRM,LTR/Ty3,LTR/Gypsy?,LTR/Gypsy-Cigr,LTR/Gypsy_Cigr,LTR/Gypsy-Troyka,LTR/Gypsy_Troyka
Bel_Pao_LTR_retrotransposon	SO:0002266	Bel_Pao_LTR_retrotransposon,LTR/Bel-Pao,LTR/RLB,Bel-Pao,Bel/Pao,LTR/BEL,LTR/Pao,LTR/Bel_Pao
Endogenous_Retrovirus_LTR_retrotransposon	SO:0002268	Endogenous_Retrovirus_LTR_retrotransposon,LTR/HERV,HERV,LTR/ERV,LTR/RLE,RLE,LTR/ERV1,LTR/ERV2,LTR/ERV3,LTR/ERV4,LTR/ERV-Foamy,LTR/ERVK,LTR/ERVL,LTR/ERV-MaLR,LTR/ERV_Foamy,LTR/ERVL_MaLR,LTR/ERVL-MaLR,Endogenous_Retrovirus,LTR/Foamy,LTR/ERV-Lenti,LTR/ERV_Lenti,LTR/Lenti
RR_tract	SO:0000435	poly_purine_tract,RR_tract
primer_binding_site	SO:0005850	primer_binding_site,PBS
long_terminal_repeat	SO:0000286	long_terminal_repeat

## nonLTR retrotransposons
non_LTR_retrotransposon	SO:0000189	non_LTR_retrotransposon,non_LTR,nonLTR/unknown,nonLTR/Unknown
LINE_element	SO:0000194	LINE_element,LINE/unknown,LINE,LINE/Unknown,LINE/NA
R2_LINE_retrotransposon	SO:0002269	R2_LINE_retrotransposon,LINE/R2,LINE/RIR,nonLTR/RIR,RIR,LINE/R2-Hero,LINE/R2-NeSL,LINE/R2_Hero,LINE/R2_NeSL,LINE/R2-Dualen,LINE/R2_Dualen
Jockey_LINE_retrotransposon	SO:0002271	Jockey_LINE_retrotransposon,LINE/Jockey,LINE/RIJ,nonLTR/RIJ,RIJ,LINE/I-Jockey
L1_LINE_retrotransposon	SO:0002272	L1_LINE_retrotransposon,LINE/L1,LINE/RIL,nonLTR/RIL,RIL,LINE-1,LINE/L1-Tx1,LINE/L1_Tx1,LINE/L1?
I_LINE_retrotransposon	SO:0002273	I_LINE_retrotransposon,LINE/I,LINE/RII,nonLTR/RII,LINE/I-Nimb,LINE/I_Nimb
RTE_LINE_retrotransposon	SO:0002270	RTE_LINE_retrotransposon,LINE/RTE,LINE/RIT,nonLTR/RIT,RIT,LINE/RTEX,LINE/RTE-X,LINE/RTE-BovB,LINE/RTE_BovB,LINE/RTE_X,LINE/RTE-RTE,LINE/RTE_RTE
SINE_element	SO:0000206	SINE_element,SINE/unknown,SINE,SINE/Unknown,SINE?,SINE?/NA,SINE/U,SINE/NA,SINE/U-L1
tRNA_SINE_retrotransposon	SO:0002274	tRNA_SINE_retrotransposon,SINE/tRNA,SINE/RST,nonLTR/RST,RST,tRNA,SINE2/tRNA,SINE/tRNA-Core-RTE,SINE/tRNA-V-CR1,tRNA,SINE/tRNA_CR1,SINE/tRNA_RTE,SINE/tRNA_V_RTE,SINE/tRNA-V-RTE,SINE/tRNA-CR1,SINE/tRNA_Core_RTE,SINE/tRNA_V_CR1,SINE/tRNA-RTE,tRNA/NA,SINE/tRNA-Deu-L2,SINE/tRNA-Mermaid,SINE/tRNA-V,SINE/tRNA-Core,SINE/tRNA-Deu,SINE/tRNA-Deu-RTE,SINE/tRNA-L1
5S_SINE_retrotransposon	SO:0002276	5S_SINE_retrotransposon,SINE/5S,SINE/RSS,nonLTR/RSS,RSS,SINE3/5S,SINE/5S-Deu-L2,SINE/5S_Deu_L2
7SL_SINE_retrotransposon	SO:0002275	7SL_SINE_retrotransposon,SINE/7SL,SINE/RSL,nonLTR/RSL,RSL,SINE1/7SL
YR_retrotransposon	SO:0002286	YR_retrotransposon,YR/unknown,YR/Unknown
Ngaro_YR_retrotransposon	SO:0002288	Ngaro_YR_retrotransposon,YR/Ngaro,YR/RYN,Ngaro,RYN,LTR/Ngaro,DIRS/Ngaro
DIRS_YR_retrotransposon	SO:0002287	DIRS_YR_retrotransposonYR/DIRS,YR/RYD,DIRS,RYD,LTR/DIRS,DIRS/NA,LTR/DIRS?,DIRS/PAT-like
Viper_YR_retrotransposon	SO:0002289	Viper_YR_retrotransposon,YR/Viper,YR/RYV,Viper,RYV
Penelope_retrotransposon	SO:0002290	Penelope_retrotransposon,Penelope,nonLTR/RPP,RPP,nonLTR/Penelope,Penelope/NA,LINE/Penelope,PLE/Chlamys

## parts
target_site_duplication	SO:0000434	target_site_duplication,TSD
U_box	SO:0001788	U_box
# END.
"""

import re


def parse_alias_file(path: str):
    """
    Returns a list of alias entries. Each entry is a dict:
      {
        'term': str,
        'so_id': str,
        'aliases': List[str],
        'first_with_slash': Optional[str],
        'leftmost': str,
        'first_with_slash_super': Optional[str],  # part after '/'
      }
    Comment lines (starting with '#') are ignored.
    """
    entries = []
    with open(path, encoding='utf-8') as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            term, so_id, alias_csv = parts[0], parts[1], parts[2]
            aliases = [a.strip() for a in alias_csv.split(',') if a.strip()]
            if not aliases:
                continue
            first_with_slash = None
            first_with_slash_super = None
            for a in aliases:
                if '/' in a:
                    first_with_slash = a
                    first_with_slash_super = a.split('/', 1)[1]
                    break
            entries.append({
                'term': term,
                'so_id': so_id,
                'aliases': aliases,
                'first_with_slash': first_with_slash,
                'leftmost': aliases[0],
                'first_with_slash_super': first_with_slash_super
            })
    return entries

# non-letter boundary (used for cautious 2–3 char matching)
_NONLETTER = r'[^A-Za-z]'

def boundary_regex(s: str) -> re.Pattern:
    """
    Build a regex that matches s only when surrounded by non-letters (or ends).
    Case-insensitive.
    """
    esc = re.escape(s)
    return re.compile(rf'(?<![A-Za-z]){esc}(?![A-Za-z])', re.IGNORECASE)

def contains_ci(hay: str, needle: str) -> bool:
    """Case-insensitive substring test."""
    return needle.lower() in hay.lower()

def eval_match(te_type: str, entry: dict) -> tuple[int, bool, bool]:
    """
    Evaluate how well te_type matches this alias entry (case-insensitive).

    Returns (score, boundary_hit, exact_super):
      - score >= 1 means a match was found
      - +2 if boundary-separated match seen
      - +3 if exact match to superfamily part of a 'subclass/superfamily' alias
    """
    aliases: list[str] = entry['aliases']
    found_any = False
    boundary_hit = False
    exact_super = False

    # exact match to superfamily (right side of first slash), case-insensitive
    sf = entry['first_with_slash_super']
    if sf and te_type.lower() == sf.lower():
        found_any = True
        exact_super = True

    # cautious for length 2–3
    if len(te_type) in (2, 3):
        bpat = boundary_regex(te_type)
        for a in aliases:
            if bpat.search(a):
                found_any = True
                boundary_hit = True
    else:
        # len >= 4: allow any substring, but still prefer boundary hits if present
        bpat = boundary_regex(te_type)
        for a in aliases:
            if contains_ci(a, te_type):
                found_any = True
                if bpat.search(a):
                    boundary_hit = True

    score = 0
    if found_any:
        score = 1
        if boundary_hit:
            score += 2
        if exact_super:
            score += 3
    return score, boundary_hit, exact_super

def choose_classification(entry: dict) -> tuple[str, bool]:
    """
    Decide classification string and whether it has a slash.
    Preference: leftmost alias containing '/', else leftmost alias (no '/').
    Returns (classification, has_slash)
    """
    if entry['first_with_slash']:
        return entry['first_with_slash'], True
    return entry['leftmost'], False

def convert_headers(
    fasta_in: str,
    so_alias_file: str,
    fasta_out: str,
    log_out: str,
    unmatched_to_unknown: bool = False
):
    entries = parse_alias_file(so_alias_file)

    n_total = 0
    n_skipped_short = 0
    n_renamed = 0
    n_unmatched = 0
    n_unmatched_fallback = 0

    with open(fasta_in, encoding='utf-8') as fin, \
         open(fasta_out, 'w', encoding='utf-8') as fout, \
         open(log_out, 'w', encoding='utf-8') as flog:

        flog.write("# Log for FASTA header conversion to RepeatMasker format\n")
        flog.write(f"# Input FASTA: {fasta_in}\n")
        flog.write(f"# Alias file:  {so_alias_file}\n")
        flog.write(f"# Option unmatched_to_unknown: {unmatched_to_unknown}\n\n")

        for line in fin:
            if not line.startswith('>'):
                fout.write(line)
                continue

            n_total += 1
            raw_header = line.rstrip('\n')[1:]  # drop '>'
            parts = raw_header.split('\t')
            if len(parts) < 2:
                # Not the expected format; leave untouched
                fout.write(line)
                flog.write(f"UNTOUCHED (format): >{raw_header}\n")
                continue

            name = parts[0].strip()
            te_type = parts[1].strip() if len(parts) >= 2 else ""
            # species not used for matching
            # species = parts[2].strip() if len(parts) >= 3 else ""

            # Skip if < 2 chars
            if len(te_type) < 2:
                n_skipped_short += 1
                fout.write(line)
                flog.write(f"SKIP short type (<2): >{raw_header}\n")
                continue

            # Find best alias match (case-insensitive logic inside eval_match)
            best = None
            best_score = 0
            best_b = False
            best_exact_sf = False

            for e in entries:
                score, b_hit, ex_sf = eval_match(te_type, e)
                if score > best_score:
                    best = e
                    best_score = score
                    best_b = b_hit
                    best_exact_sf = ex_sf

            if not best or best_score == 0:
                # No alias hit
                n_unmatched += 1
                if unmatched_to_unknown:
                    # Force to unknown/unknown (has '/'), keep original name as ID
                    new_header = f">{name}#unknown/unknown"
                    fout.write(new_header + "\n")
                    n_unmatched_fallback += 1
                    flog.write(
                        f"FALLBACK unknown/unknown: >{raw_header}  ->  {new_header}\n"
                    )
                else:
                    # Leave untouched
                    fout.write(line)
                    flog.write(f"UNMATCHED (no alias hit): >{raw_header}\n")
                continue

            # Build classification
            classification, has_slash = choose_classification(best)

            # with '/', keep original name; without '/', adopt alias as ID
            if has_slash:
                new_id = name
                trailing = classification
            else:
                new_id = classification
                trailing = classification

            new_header = f">{new_id}#{trailing}"
            fout.write(new_header + "\n")
            n_renamed += 1

            why = []
            if best_exact_sf:
                why.append("exact_superfamily")
            if best_b:
                why.append("boundary")
            if not why:
                why.append("substring")

            flog.write(
                "RENAME: "
                f">{raw_header}  ->  {new_header}  "
                f"[term={best['term']}; score={best_score}; why={'+'.join(why)}]\n"
            )

        # Summary
        flog.write("\n# Summary\n")
        flog.write(f"Total headers:          {n_total}\n")
        flog.write(f"Renamed:                {n_renamed}\n")
        flog.write(f"Skipped (<2):           {n_skipped_short}\n")
        flog.write(f"Unmatched:              {n_unmatched}\n")
        flog.write(f"Fallback unknown/unknown applied: {n_unmatched_fallback}\n")

def main():
    import argparse
    p = argparse.ArgumentParser(
        description=(
            "Convert FASTA headers (name\\ttype\\tspecies) to RepeatMasker format "
            "using an SO alias list. Matching is case-insensitive. "
            "Short types (<2) are skipped. 2–3 char types require non-letter boundaries."
        )
    )
    p.add_argument("-i", "--in-fasta", required=True,
                   help="Input FASTA with headers 'name<TAB>type<TAB>species'")
    p.add_argument("-k", "--key", required=True,
                   help="Alias key file (SO.txt-like; columns: term, SO_ID, aliases)")
    p.add_argument("-o", "--out-fasta", required=True,
                   help="Output FASTA with RepeatMasker-style headers")
    p.add_argument("-l", "--log", required=True,
                   help="Log file for skips/renames")
    p.add_argument("--unmatched-to-unknown", action="store_true",
                   help="If set, dump UNMATCHED to '>name#unknown/unknown' instead of leaving unchanged.")
    args = p.parse_args()

    convert_headers(
        args.in_fasta,
        args.key,
        args.out_fasta,
        args.log,
        unmatched_to_unknown=args.unmatched_to_unknown
    )

if __name__ == "__main__":
    main()
