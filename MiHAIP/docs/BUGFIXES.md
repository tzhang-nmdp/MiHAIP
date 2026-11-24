# MiHAIP (MiHAIP-master) — Bug fixes in the Java source
Documents two bugs found while checking why missense variants were missing from
`protein/*_protien.fasta` and which FASTA record is WT. Both are fixed in the
`.java` files in `MiHAIP/src`.
> Source note: code snippets are taken from the current `.java` files in
> `MiHAIP/src` (`database/ImportData`, `analysis/{Transcript, ProcessVcf,
> GenerateOutput}`). "Before" snippets are reconstructed from the investigation
> reports in `MiHAIP/.agents/tasks/` (`missing-missense-in-fasta.md`,
> `sample52-merge-vs-old-fasta.md`, `wt-vs-alt-record-order.md`) and are not
> from a version-control diff (the folder is not a git repository). Numbers come
> from those reports and from the scratch scripts in `%TEMP%\mihaip_strand`
> (`verify.py`, `StrandTest.java`). **None of the fixes has been run end to end:
> no JVM run, no `cds.db` rebuild, no new FASTA.** See the Verification section
> of each bug.
---
## 0. Overview
| # | Bug | File / method | Effect before the fix | Fix | Needs |
|---|---|---|---|---|---|
| 1 | Transcript version cut to one digit | `database/ImportData` (`idPattern`); warning in `analysis/GenerateOutput.processGene` | `ENST….10` and higher stored as `.1`, `getSequence` returns `""`, variant silently missing from the FASTA | `\d` to `\d+` in the regex; stderr warning when a sequence is missing | Rebuild `cds.db`, rerun Merge |
| 2 | Minus-strand alleles written unchanged into the CDS | `analysis/Transcript` (strand), `analysis/ProcessVcf.processLine/processTrans`, `analysis/GenerateOutput.processGene` | Wrong or swapped WT and mutant residues for minus-strand genes | Detect strand from HGVS.c, complement ref/alt before writing into the CDS | Rerun Merge and everything downstream |
Design that bug 2 restores: **FASTA record 1 = WT (VCF REF), record 2 = mutant
(VCF ALT)**, for both plus- and minus-strand genes.
## 1. How to apply
1. **Port the fixes to the real source first.** `MiHAIP/src` is not the code
   that produced the deployed output.
   - The deployed jar (`Tools/versions/MiHAIP0.3.3_patched2.jar`, run by
     `Tools/run_miha.sh`) differs from `src`: `MergeTool` (same pos and same alt
     removes both records in the jar, `src` keeps one record with ref and alt
     swapped), `ProcessVcf` (jar groups consecutive lines of one transcript,
     `src` makes one `Transcript` per line), and `GenerateOutput` (`patched2` has
     an older build without the 41-aa window; the sample FASTAs match
     `MiHAIP0.3.3_new.jar`).
   - `src` now compiles (0 errors, jdk21 javac, classpath
     `sqlite-jdbc-3.8.11.2.jar`). `analysis/Patient.java` was missing from every
     tree and jar, so it was reconstructed from usage (enum `Donor`/`Reciepor`
     with `getID()` returning "D"/"R"; the real ID values are unknown). The
     stray `analysis/Transcript_src.java` was moved to
     `../backup/Transcript_src.java.bak`.
   - So the changes below must be re-applied to the source tree the jar is really
     built from, then the jar rebuilt.
2. **Compile and package** the jar from that tree (for example with the bundled
   `jdk21` under `MiHAIP_1.4.5`).
3. **Rebuild `cds.db`** (bug 1): delete the old `Tools/cds.db` and rerun the CDS
   import from `Homo_sapiens.GRCh38.cds.all.fa` with the fixed `ImportData`.
   That FASTA was not on the machine when this was checked, so the rebuild has
   not been done.
4. **Rerun Step 4 (Merge)** for every sample, then all downstream steps (netChop,
   chop, NetMHCpan, affinity/frequency, compare). Old `protein/`, `DNA/` and
   `meta/` files are stale.
5. Look at stderr for `WARNING: no CDS sequence found for transcript …`. Each
   line is a variant that still cannot be placed (see section 4).
---
## 2. Bug 1 — Transcript version truncated in `cds.db`
### 2.1 Summary
`ImportData` stored transcript IDs with only one digit after the dot. SnpEff
reports `ENST00000258381.10`; the database held `ENST00000258381.1`. The exact
match in `DatabaseUtil.getSequence` failed, and `GenerateOutput.processGene`
wrote nothing for that transcript.
### 2.2 Symptom with example
Sample 1: 20 missense VCF records (17 in `1_fp_ann_msv.vcf`, 3 in
`1_fn_ann_msv.vcf`), 17 reached `1_protien.fasta`. The 3 lost:
| Variant | Source | Gene | In VCF | In `cds.db` |
|---|---|---|---|---|
| chr2:230207994 C>G p.Gly299Arg | fp | SP110 | `ENST00000258381.10` | `ENST00000258381.1` |
| chr6:143502020 T>C p.Met356Val | fp | FUCA2 | `ENST00000002165.10` | `ENST00000002165.1` |
| chr7:44977022 C>T p.Val49Met | fp | MYO1G | `ENST00000258787.11` | `ENST00000258787.1` |
Sample 52: 7,584 ref/alt pairs in `52_protien_old.fasta`. After the fix the
emulated Merge gives 8,274 pairs: **690 pairs from 420 transcripts** were lost
(versions `.10` to `.16`: `.10` 383 pairs, `.11` 181, `.12` 79, `.13` 42, `.15` 4,
`.16` 1). Nothing in the old FASTA was removed or changed by the fix.
### 2.3 Root cause
`database/ImportData.java`, field `idPattern`: the regex `\d` matched one digit
after the dot. `DatabaseUtil.getSequence` runs `WHERE TRANSCRIPT_ID = '<id>'`
(exact match), returns `""` when no row is found, and nothing is logged. In
`processGene` the loop `if(i < dna.length())` is then never true, so no record
is written.
### 2.4 Fix
`database/ImportData.java`:
```java
// before
static String idPattern = ">ENST\\d+\\.\\d";
// after
static String idPattern = ">ENST\\d+\\.\\d+";
```
`analysis/GenerateOutput.java`, start of `processGene` (new warning):
```java
String dna = DatabaseUtil.getSequence(gene.transcriptID);
if(dna == null || dna.isEmpty()){
    System.err.println("WARNING: no CDS sequence found for transcript " + gene.transcriptID + ", variants skipped");
}
```
The warning does not change control flow. An empty sequence still produces
no FASTA record, as before.
### 2.5 Verification
Checked (emulation in Python, no JVM, `cds.db` not rebuilt):
- The emulation of the deployed Merge logic on the unchanged `cds.db` reproduced
  `52_protien_old.fasta` exactly (7,584 pairs, same headers, sequences, order).
- "Fixed" mode takes the sequence from the truncated row (`….N` for `….NN`) and
  accepts it only if its length equals the CDS length in ANN field 13. All 420
  recovered transcripts pass. The CDS base at the HGVS.c position matches the VCF
  ref (or its complement) for 690 of 691 recovered variants.
- Old FASTA is an ordered subsequence of the emulated new output (old-only = 0).
Not verified:
- No real `cds.db` rebuild, no run of the Java code, no check that the new
  warning prints.
- A real rebuild may differ if the Ensembl CDS release differs from the SnpEff
  release (this is the MPRIP case, section 5).
### 2.6 Impact on outputs
More records in `protein/<id>_protien.fasta`, `DNA/*`, `meta/*` and everything
downstream, for every transcript with version 10 or higher. Existing records do
not change. New records of the recovered set are still subject to bug 2 until it
is applied.
### 2.7 Required actions
Apply to the real source, rebuild `cds.db`, rerun Merge and downstream (section 1).
---
## 3. Bug 2 — Strand handling in `GenerateOutput`
### 3.1 Summary
VCF REF/ALT are bases of the forward genomic strand. The CDS from `cds.db` is in
transcript orientation. `processGene` wrote the VCF bases into the CDS without
reverse complementing for minus-strand transcripts, so the substituted base was
wrong. Record 1 (should be WT) and record 2 (should be mutant) got wrong
residues.
### 3.2 Symptom with example
Sample 1, `1_fp_ann_msv.vcf` (GT `1/1`): chr5:96803547 C>G, ERAP1,
`ENST00000296754.7`, `c.380G>C p.Arg127Pro`. HGVS.c ref `G` is the complement of
VCF ref `C`, so the gene is on the minus strand.
```
Before the fix
>chr5|96803547|ERAP1|ENST00000296754.7|127
LRKGAGERLSEEPLQVLEHPpQEQIALLAPEPLLVGLPYTV      record 1 = Pro (mutant)
>chr5|96803547|ERAP1|ENST00000296754.7|127
LRKGAGERLSEEPLQVLEHPrQEQIALLAPEPLLVGLPYTV      record 2 = Arg (WT)
```
Expected after the fix: Arg on record 1, Pro on record 2.
Plus-strand control, chr6:12123016 A>G, HIVEP1, `c.3221A>G p.Asn1074Ser`:
record 1 `n`, record 2 `s`. Correct before and after.
Why ERAP1 was an exact swap: for minus-strand `C>G`, `G>C`, `A>T`, `T>A`, the
genomic ALT equals the complement of REF, which is the transcript WT base, and
vice versa. For minus-strand transitions (`C>T`, `G>A`, `A>G`, `T>C`) the written
base matches neither allele, so both records are wrong (TYMP `c.1408G>A
p.Ala470Thr` gave `p`/`s`; CRLF2 `c.968C>T p.Ser323Phe` gave `c`/`y`).
Note that the pair is not a heterozygous-genotype effect: the code always builds
one ref sequence and one alt sequence from one VCF record and never reads GT.
### 3.3 Root cause
`analysis/GenerateOutput.processGene` replaced the CDS base with `gene.ref` (for
record 1, since `Merge` calls it with `replaceTwo = true`) and `gene.alt` (record
2) as they were in the VCF. No reverse complement existed anywhere in the code.
`gene.ref` and `gene.alt` are VCF columns 4 and 5, set by
`ProcessVcf.processLine` via `gene.setNT(basic[3], basic[4])`.
### 3.4 Fix
**How strand is determined.** Not from the database. `Transcript.setStrandFromHgvsC`
reads the SnpEff HGVS.c field (ANN field 10, e.g. `c.380G>C`) and compares its
reference base with the VCF REF base:
- equal: plus strand
- equal to the complement of the VCF REF: minus strand
- anything else (indels, no single-base substitution, REF longer than one base,
  null): the current value is kept (default plus)
`analysis/Transcript.java` (new field, copy in `setInfo`, new methods):
```java
//true if the transcript is on the minus strand (CDS is the reverse complement of the VCF bases)
boolean minusStrand = false;

public void setInfo(Transcript t){
    ...
    this.minusStrand = t.minusStrand;   // added
}

public void setStrandFromHgvsC(String hgvsC){
    if(hgvsC == null || ref == null || ref.length() != 1){ return; }
    Matcher m = Pattern.compile("c\\.\\d+([ACGT])>[ACGT]").matcher(hgvsC);
    if(!m.find()){ return; }
    String cRef = m.group(1);
    String vcfRef = ref.toUpperCase();
    if(cRef.equals(vcfRef)){ minusStrand = false; }
    else if(cRef.equals(complement(vcfRef))){ minusStrand = true; }
}

public static String complement(String s){ ... }   // A<->T, C<->G, case kept, other chars unchanged
```
`analysis/ProcessVcf.java` (calls added so that every `Transcript` gets a strand;
`setStrandFromHgvsC` needs `ref`, so it runs after `setNT`):
```java
// processLine, after gene.addPos(DNAchange, gene.getAlt());
gene.setStrandFromHgvsC(detail[9]);

// processTrans, for the additional transcript created from the same record
Transcript newGene = new Transcript(transcript);
newGene.setInfo(gene);                       // copies minusStrand too
newGene.addPos(DNAchange, gene.getAlt());
newGene.setStrandFromHgvsC(detail[3]);       // added: strand of this transcript
geneList.add(newGene);
```
`analysis/GenerateOutput.java`, `processGene`:
```java
// before
sbRef.replace(i-1, i, gene.ref.toLowerCase());     // replaceTwo
sbAlt.replace(i-1, i, gene.alt.toLowerCase());

// after
String refNT = gene.ref;
String altNT = gene.alt;
if(gene.minusStrand){
    refNT = Transcript.complement(refNT);
    altNT = Transcript.complement(altNT);
}
...
sbRef.replace(i-1, i, refNT.toLowerCase());        // when replaceTwo, else CDS base as before
sbAlt.replace(i-1, i, altNT.toLowerCase());
```
Metadata file (`meta/<id>_meta.txt`). It is used for the frequency lookup by
genomic base, so it must stay on the VCF strand. The bases read back from the
CDS strings are complemented again for minus-strand genes:
```java
String metaRef = sbRef.toString().substring(DNAchangePos.get(i)-1, DNAchangePos.get(i));
String metaAlt = sbAlt.toString().substring(DNAchangePos.get(i)-1, DNAchangePos.get(i));
if(gene.minusStrand){
    metaRef = Transcript.complement(metaRef);
    metaAlt = Transcript.complement(metaAlt);
}
```
The refDNA and altDNA dumps now hold the transcript-sense bases.
### 3.5 Verification
Checked:
- **Strand rule** (`%TEMP%\mihaip_strand\StrandTest.java`, compiled class
  present): written to test `C`/`c.380G>C`, `A`/`c.3221A>G`, `G`/`c.1408G>A`,
  `G`/`c.1013C>T`, `A`/`c.5_6del` and `complement("AcGt")`. The compiled class is
  there; the printed output is not saved in the folder, so this document does not
  claim a result.
- **Emulation on the old output** (`verify.py`, Python; reads the FASTA, the VCFs
  and the CDS codon tables, no Java). Sample 1: old output on the 15 uniquely
  matched variants: 8 minus-strand inconsistent, 7 plus-strand consistent, 2 not
  matched; emulated new logic: 7 minus-strand consistent, 7 plus-strand
  consistent, 1 minus-strand emulation failed.
- **Sample 52** (matched single-base missense variants whose old output could be
  checked against HGVS.p): old output 1,978 of 5,137 pairs consistent (about
  38%; 1,901 plus, 77 minus), and 2,522 minus-strand pairs inconsistent. Emulated
  new logic: 1,273 minus-strand consistent, 616 ambiguous (several
  codon-context candidates, some consistent), 9 inconsistent, 701 not
  reproducible from the old output. Plus-strand is unchanged (1,901 consistent,
  637 inconsistent; these 637 come from grouped transcripts and multi-allelic
  records in the old jar logic, not from strand).
- The 25 variants whose HGVS.c base matches neither VCF REF nor its complement
  (sample 52 old output) keep the default plus strand.
Not verified:
- **The Java was compiled but never run** (compile check only, section 1).
  Everything above is emulation. `StrandTest` only covers `Transcript`.
- No new FASTA was generated, and the real result for the 701 non-reproducible
  and 616 ambiguous pairs is unknown until Merge is rerun.
- The fixed code is not in the deployed jar. The jar's grouped `Transcript`
  objects reuse the last line's alleles (section 5), so the same rerun still
  gives wrong residues for grouped transcripts.
### 3.6 Impact on outputs
Sample 52 old output: 3,732 variants plus strand, 3,827 minus strand, 25 neither
(CDS base at the HGVS.c position matches VCF ref, its complement, or neither).
Only about 38% of pairs had correct ref/alt residues before the fix.
- `protein/*.fasta`: residues change for minus-strand genes. Record 1 is WT and
  record 2 is mutant for both strands.
- Peptides from NetChop, NetMHCpan and compare: peptide to allele assignment
  changes for minus-strand genes. No peptide is dropped by the swap itself, but
  peptide sets change, and results from the old runs are not valid for these
  genes. Nothing downstream labels the sets WT or mutant; they are treated
  symmetrically.
- `DNA/*_refDNA.txt` and `*_altDNA.txt`: minus-strand lines are transcript-sense.
- `meta/*`: the same genomic REF/ALT bases as before (kept on the VCF strand).
### 3.7 Required actions
Apply to the real source, rerun Merge and all downstream steps (section 1). Spot
check ERAP1 in sample 1: record 1 should show `r`, record 2 `p`.
---
## 4. Messages and detection after the fixes
- `WARNING: no CDS sequence found for transcript <id>, variants skipped` on stderr
  marks every transcript that still has no CDS row.
- The old FASTA files are unchanged. Compare new and old with the header
  (`>chr|pos|GENE|ENST|aa`); X is written as `chr24` and Y as `chr25`.
---
## 5. Remaining known issues (not fixed)
| Issue | Count / example | Why it is open |
|---|---|---|
| First ANN is `structural_interaction_variant` or `protein_protein_contact`, whose transcript field is a PDB ID such as `3T6A:A_598-A_645:ENST00000260502` | 65 variants in sample 52 skipped | `ProcessVcf.processLine` uses `detail[6]` of the first ANN only; `getSequence` returns `""`. All 65 have a later `missense_variant` ANN whose transcript is in `cds.db` (for example BCAR3 chr1:93571710 T>C, `ENST00000260502.10`). Fix would go in `ProcessVcf`. |
| Transcript version mismatch: MPRIP chr17:17165489 A>G p.Arg668Gly, VCF `ENST00000313485.10` | 1 variant in sample 52 | `cds.db` has `ENST00000313485.1` with 5,385 nt, SnpEff CDS length is 5,382, so it is a different transcript version. Not fixed by bug 1 unless the CDS FASTA release matches SnpEff. |
| Secondary ANN entries are ignored | 212 records in sample 52 have more than one missense ANN | `ProcessVcf.processTrans` never matches (later entries begin with the allele, not `ENST`), so secondary transcripts do not reach the FASTA. |
| `MergeTool` genotype-mismatch behaviour differs from `STEP1.md` | 2,163 FP and 2,163 FN records cancelled in sample 52 (4,326 lines) | The deployed jar removes a record present with the same position and alt in both files. `STEP1.md` describes keeping the donor record. `src/MergeTool` keeps one record (the FN one) with ref and alt swapped; the jar does not. The docs and the code disagree; no decision was made on which is intended. |
| Grouped transcripts (jar only) | Sample 1: GTPBP6 `m`/`t` at Thr338Met; sample 52: grouped records | The jar's `ProcessVcf` merges consecutive lines of one transcript and writes the last line's ref/alt at every position, and `aachangePos` follows `HashMap` order so `Translator` may not lowercase some residues. `src` makes one `Transcript` per line, so it does not have this bug. Still open in the deployed jar, and its strand fix has the same limit. |
| ref window equals alt window | 297 old pairs in sample 52 | Seen in the old output, in part a consequence of the strand bug and the grouped alleles. Not rechecked after the fix. |
| Stop codons in output | 555 old pairs in sample 52 have `*`; 6 have `X` | Not investigated. |

---
## 6. Port to `Tools/versions/MiHAIP0.3.3_new` (source tree)
The fixes were also applied to the decompiled tree `Tools/versions/MiHAIP0.3.3_new`.
- **Changed:** `analysis/Transcript.java` only (`minusStrand`, copy in `setInfo`,
  `setStrandFromHgvsC`, static `complement`, plus the `Matcher`/`Pattern` imports).
  Style of the decompiled file kept. No `Patient` field was added, so
  `analysis/Patient.java` is not needed in that tree (the earlier `cannot find
  symbol Patient` error came from compiling the MiHAIP-master `Transcript.java`
  inside it).
- **Already fixed in that tree before this port (no edit made):**
  `analysis/ProcessVcf.java` (`setStrandFromHgvsC` calls),
  `analysis/GenerateOutput.java` (complement of ref/alt, meta kept on VCF strand,
  missing-CDS warning) and `database/ImportData.java` (`\d+` in `idPattern`).
  They did not compile without the new `Transcript` members.
- **Backup:** `MiHAIP0.3.3_new/backup_prefix/Transcript.java` (pre-port file).
- **Compile:** jdk21 `javac -encoding UTF-8 -nowarn -cp
  "database\sqlite-jdbc-3.15.0.jar;MiHAIP0.3.3_new.jar" -d %TEMP%\mihaip033_build
  analysis\Transcript.java analysis\ProcessVcf.java analysis\GenerateOutput.java
  analysis\ProcessAff.java database\ImportData.java`: 0 errors. Other classes
  come from the jar (the tree has only a few `.java` files; `clevage/Clevage.java`
  has a BOM and was not compiled).
- **Run check:** a throwaway harness (`ProcessVcf` + `GenerateOutput`, real
  `Tools/cds.db`, sample-1 `1_fp/1_fn_ann_msv.vcf`, output in `%TEMP%`) gave ERAP1
  `...HP`**`r`**`QEQ...` on record 1 and `p` on record 2, HIVEP1 `n` then `s`. The
  3 sample-1 transcripts with version 10+ printed the missing-CDS warning because
  `cds.db` is still the old truncated one.
- **Still to do:** rebuild the jar from this tree, rebuild `cds.db`, rerun Merge.
