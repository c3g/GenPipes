package ca.mcgill.genome.bamtools.maf;

import java.io.Closeable;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.EnumMap;
import java.util.List;

import ca.mcgill.genome.bamtools.commands.VCF2MAF;
import ca.mcgill.genome.bamtools.parsers.GeneMapper;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect.FormatVersion;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

/**
 * https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification
 * @author lletourn
 *
 */
public class MAFWriter implements Closeable {
	private final PrintStream writer;
	private String normal;
	private String tumor;
	private GeneMapper mapper;

	private enum VariantType {
		SNP, DNP, TNP, ONP, INS, DEL, Consolidated
	};

	private static final EnumMap<EffectType, Integer> effectPriority = new EnumMap<EffectType, Integer>(EffectType.class);
	private static final EnumMap<EffectType, String> snpEff2MAF = new EnumMap<EffectType, String>(EffectType.class);
	static {
		effectPriority.put(EffectType.RARE_AMINO_ACID, 1);
		effectPriority.put(EffectType.SPLICE_SITE_ACCEPTOR, 2);
		effectPriority.put(EffectType.SPLICE_SITE_DONOR, 3);
		effectPriority.put(EffectType.START_LOST, 4);
		effectPriority.put(EffectType.EXON_DELETED, 5);
		effectPriority.put(EffectType.FRAME_SHIFT, 6);
		effectPriority.put(EffectType.STOP_GAINED, 7);
		effectPriority.put(EffectType.STOP_LOST, 8);
		effectPriority.put(EffectType.NON_SYNONYMOUS_CODING, 9);
		effectPriority.put(EffectType.CODON_CHANGE, 10);
		effectPriority.put(EffectType.CODON_INSERTION, 11);
		effectPriority.put(EffectType.CODON_CHANGE_PLUS_CODON_INSERTION, 12);
		effectPriority.put(EffectType.CODON_DELETION, 13);
		effectPriority.put(EffectType.CODON_CHANGE_PLUS_CODON_DELETION, 14);
		effectPriority.put(EffectType.UTR_5_DELETED, 15);
		effectPriority.put(EffectType.UTR_3_DELETED, 16);
		effectPriority.put(EffectType.SYNONYMOUS_START, 17);
		effectPriority.put(EffectType.NON_SYNONYMOUS_START, 18);
		effectPriority.put(EffectType.START_GAINED, 19);
		effectPriority.put(EffectType.SYNONYMOUS_CODING, 20);
		effectPriority.put(EffectType.SYNONYMOUS_STOP, 21);
		effectPriority.put(EffectType.NON_SYNONYMOUS_STOP, 22);
		effectPriority.put(EffectType.UTR_5_PRIME, 23);
		effectPriority.put(EffectType.UTR_3_PRIME, 24);
		effectPriority.put(EffectType.REGULATION, 25);
		effectPriority.put(EffectType.UPSTREAM, 26);
		effectPriority.put(EffectType.DOWNSTREAM, 27);
		effectPriority.put(EffectType.GENE, 28);
		effectPriority.put(EffectType.TRANSCRIPT, 29);
		effectPriority.put(EffectType.EXON, 30);
		effectPriority.put(EffectType.INTRON_CONSERVED, 31);
		effectPriority.put(EffectType.INTRON, 32);
		effectPriority.put(EffectType.INTRAGENIC, 33);
		effectPriority.put(EffectType.INTERGENIC, 34);
		effectPriority.put(EffectType.INTERGENIC_CONSERVED, 35);
		effectPriority.put(EffectType.NONE, 36);
		effectPriority.put(EffectType.CHROMOSOME, 37);
		effectPriority.put(EffectType.CUSTOM, 38);
		effectPriority.put(EffectType.CDS, 39);

		snpEff2MAF.put(EffectType.NONE, "");
		snpEff2MAF.put(EffectType.CHROMOSOME, "");
		snpEff2MAF.put(EffectType.CUSTOM, "");
		snpEff2MAF.put(EffectType.GENOME, "");
		snpEff2MAF.put(EffectType.MICRO_RNA, "");
		snpEff2MAF.put(EffectType.INTERGENIC, "IGR");
		snpEff2MAF.put(EffectType.UPSTREAM, "5'Flank");
		snpEff2MAF.put(EffectType.UTR_5_PRIME, "5'UTR");
		snpEff2MAF.put(EffectType.UTR_5_DELETED, "5'UTR");
		snpEff2MAF.put(EffectType.START_GAINED, "De_novo_Start_InFrame");
		snpEff2MAF.put(EffectType.SPLICE_SITE_ACCEPTOR, "Splice_Site");
		snpEff2MAF.put(EffectType.SPLICE_SITE_BRANCH, "");
		snpEff2MAF.put(EffectType.SPLICE_SITE_BRANCH_U12, "");
		snpEff2MAF.put(EffectType.SPLICE_SITE_DONOR, "Splice_Site");
		snpEff2MAF.put(EffectType.START_LOST, "Missense_Mutation");
		snpEff2MAF.put(EffectType.SYNONYMOUS_START, "Silent");
		snpEff2MAF.put(EffectType.NON_SYNONYMOUS_START, "Missense_Mutation");
		snpEff2MAF.put(EffectType.CDS, "Targeted_Region");
		snpEff2MAF.put(EffectType.GENE, "Targeted_Region");
		snpEff2MAF.put(EffectType.TRANSCRIPT, "RNA");
		snpEff2MAF.put(EffectType.EXON, "Targeted_Region");
		snpEff2MAF.put(EffectType.EXON_DELETED, "Frame_Shift_Del");
		snpEff2MAF.put(EffectType.NON_SYNONYMOUS_CODING, "Missense_Mutation");
		snpEff2MAF.put(EffectType.SYNONYMOUS_CODING, "Silent");
		snpEff2MAF.put(EffectType.CODON_CHANGE, "Missense_Mutation");
		snpEff2MAF.put(EffectType.CODON_INSERTION, "In_Frame_Ins");
		snpEff2MAF.put(EffectType.CODON_CHANGE_PLUS_CODON_INSERTION, "In_Frame_Ins");
		snpEff2MAF.put(EffectType.CODON_DELETION, "In_Frame_Del");
		snpEff2MAF.put(EffectType.CODON_CHANGE_PLUS_CODON_DELETION, "In_Frame_Del");
		snpEff2MAF.put(EffectType.RARE_AMINO_ACID, "Missense_Mutation");
		snpEff2MAF.put(EffectType.STOP_GAINED, "Nonsense_Mutation");
		snpEff2MAF.put(EffectType.SYNONYMOUS_STOP, "Silent");
		snpEff2MAF.put(EffectType.NON_SYNONYMOUS_STOP, "Nonsense_Mutation");
		snpEff2MAF.put(EffectType.STOP_LOST, "Nonstop_Mutation");
		snpEff2MAF.put(EffectType.INTRON, "Intron");
		snpEff2MAF.put(EffectType.UTR_3_PRIME, "3'UTR");
		snpEff2MAF.put(EffectType.UTR_3_DELETED, "3'UTR");
		snpEff2MAF.put(EffectType.DOWNSTREAM, "3'Flank");
		snpEff2MAF.put(EffectType.INTRON_CONSERVED, "Intron");
		snpEff2MAF.put(EffectType.INTERGENIC_CONSERVED, "IGR");
		snpEff2MAF.put(EffectType.INTRAGENIC, "Targeted_Region");
		snpEff2MAF.put(EffectType.REGULATION, "5'Flank");
		snpEff2MAF.put(EffectType.FRAME_SHIFT, "Frame_Shift_");//Frame_Shift_Ins, Frame_Shift_Del
	}

	public MAFWriter(PrintStream writer) {
		this.writer = writer;
	}

	public void write(VcfEntry vcfEntry) {
		// MAF 2.3, only somatic calls are valid.

		VariantType variantType = VariantType.SNP;
		int refLength = vcfEntry.getRef().length();
		int alt1Length = vcfEntry.getAlts()[0].length();
		if (refLength > 1 || alt1Length > 1) {
			if (refLength > alt1Length) {
				variantType = VariantType.DEL;
			} else if (refLength < alt1Length) {
				variantType = VariantType.INS;
			} else {
				if (refLength == 2) {
					variantType = VariantType.DNP;
				} else if (refLength == 3) {
					variantType = VariantType.TNP;
				} else {
					variantType = VariantType.ONP;
				}
			}
		}

		String variantid = vcfEntry.getId();
		if (variantid == null || !variantid.startsWith("rs")) {
			variantid = "novel";
		}

		List<VcfEffect> allEffects = vcfEntry.parseEffects(FormatVersion.FORMAT_SNPEFF_3);
		List<VcfEffect> effects = new ArrayList<VcfEffect>();

		if(mapper == null){
			effects = allEffects;
		}
		else {
			for(VcfEffect effect : allEffects) {
				if(effect.getTranscriptId() == null) {
					effects.add(effect);
				}
				else if(mapper.getTranscriptId(effect.getTranscriptId()) != null) {
					effects.add(effect);
				}
			}

			if(effects.size() == 0) {
				VcfEffect effect = new VcfEffect("INTERGENIC(MODIFIER||||||||)", FormatVersion.FORMAT_SNPEFF_3);
				effects.add(effect);
			}
		}

		if (effects.size() > 0) {
			Collections.sort(effects, new Comparator<VcfEffect>() {
				public int compare(VcfEffect o1, VcfEffect o2) {
					return effectPriority.get(o1.getEffect()).compareTo(effectPriority.get(o2.getEffect()));
				}
			});
			VcfEffect effect = effects.get(0);
//			for (VcfEffect effect : effects) {
			String effectType = snpEff2MAF.get(effect.getEffect());

			if (effectType.equals("Frame_Shift_")) {
				if (refLength < alt1Length) effectType = "Frame_Shift_Ins";
				else effectType = "Frame_Shift_Del";
			}

			if (mapper != null) {
				if(effect.getTranscriptId() != null && mapper.getEntrezGeneId(effect.getTranscriptId()) == null) {
					System.err.println("Entrez not found for: " + effect.getTranscriptId());
				}
				writer.print(mapper.getHUGONomenclatureSymbol(effect.getTranscriptId())); //"Hugo_Symbol"
				writer.print('\t');
				writer.print(mapper.getEntrezGeneId(effect.getTranscriptId()));//"Entrez_Gene_Id"
				writer.print('\t');
			} else {
				writer.print('\t');
				writer.print('\t');
			}
			writer.print("CNG"); // Center
			writer.print('\t');
			writer.print("GRCh37"); //NCBI_Build
			writer.print('\t');
			writer.print(vcfEntry.getChromosomeName());// "Chromosome"
			writer.print('\t');
			writer.print(vcfEntry.getStart()+1);// "Start_Position");
			writer.print('\t');
			writer.print(vcfEntry.getStart()+1 + refLength - 1);//"End_Position");
			writer.print('\t');
			writer.print('+');//"Strand");
			writer.print('\t');
			writer.print(effectType);//"Variant_Classification");
			writer.print('\t');
			writer.print(variantType);//"Variant_Type");
			writer.print('\t');
			writer.print(vcfEntry.getRef());//"Reference_Allele");
			writer.print('\t');
			String genotype = vcfEntry.getVcfGenotype(VCF2MAF.TUMOR_VCF_IDX).getGenotypeStr();
			writer.print(genotype.substring(0, genotype.indexOf("/")));//"Tumor_Seq_Allele1");
			writer.print('\t');
			writer.print(genotype.substring(genotype.indexOf("/")+1));//"Tumor_Seq_Allele2");
			writer.print('\t');
			writer.print(variantid);//"dbSNP_RS");
			writer.print('\t');
			writer.print("");//"dbSNP_Val_Status");
			writer.print('\t');
			writer.print(tumor);//"Tumor_Sample_Barcode");
			writer.print('\t');
			writer.print(normal);//"Matched_Norm_Sample_Barcode");
			writer.print('\t');
			genotype = vcfEntry.getVcfGenotype(VCF2MAF.NORMAL_VCF_IDX).getGenotypeStr();
			writer.print(genotype.substring(0, genotype.indexOf("/")));//"Match_Norm_Seq_Allele1");
			writer.print('\t');
			writer.print(genotype.substring(genotype.indexOf("/")+1));//"Match_Norm_Seq_Allele2");
			writer.print('\t');
			writer.print("");//"Tumor_Validation_Allele1");
			writer.print('\t');
			writer.print("");//"Tumor_Validation_Allele2");
			writer.print('\t');
			writer.print("");//"Match_Norm_Validation_Allele1");
			writer.print('\t');
			writer.print("");//"Match_Norm_Validation_Allele2");
			writer.print('\t');
			writer.print("");//"Verification_Status");
			writer.print('\t');
			writer.print("");//"Validation_Status");
			writer.print('\t');
			writer.print("Somatic");//"Mutation_Status");
			writer.print('\t');
			writer.print("");//"Sequencing_Phase");
			writer.print('\t');
			writer.print("");//"Sequence_Source");
			writer.print('\t');
			writer.print("");//"Validation_Method");
			writer.print('\t');
			writer.print("");//"Score");
			writer.print('\t');
			writer.print("");//"BAM_File");
			writer.print('\t');
			writer.print("Illumina HiSeq");//"Sequencer");
			writer.print('\t');
			writer.print("");//"Tumor_Sample_UUID");
			writer.print('\t');
			writer.print("");//"Matched_Norm_Sample_UUID");
			writer.println();
//			}
		}
	}

	public void close() throws IOException {
		writer.close();
	}

	public GeneMapper getMapper() {
		return mapper;
	}

	public void setMapper(GeneMapper mapper) {
		this.mapper = mapper;
	}

	public String getNormal() {
		return normal;
	}

	public void setNormal(String normal) {
		this.normal = normal;
	}

	public String getTumor() {
		return tumor;
	}

	public void setTumor(String tumor) {
		this.tumor = tumor;
	}

	public void printHeader() {
		writer.print("Hugo_Symbol");
		writer.print('\t');
		writer.print("Entrez_Gene_Id");
		writer.print('\t');
		writer.print("Center");
		writer.print('\t');
		writer.print("NCBI_Build");
		writer.print('\t');
		writer.print("Chromosome");
		writer.print('\t');
		writer.print("Start_Position");
		writer.print('\t');
		writer.print("End_Position");
		writer.print('\t');
		writer.print("Strand");
		writer.print('\t');
		writer.print("Variant_Classification");
		writer.print('\t');
		writer.print("Variant_Type");
		writer.print('\t');
		writer.print("Reference_Allele");
		writer.print('\t');
		writer.print("Tumor_Seq_Allele1");
		writer.print('\t');
		writer.print("Tumor_Seq_Allele2");
		writer.print('\t');
		writer.print("dbSNP_RS");
		writer.print('\t');
		writer.print("dbSNP_Val_Status");
		writer.print('\t');
		writer.print("Tumor_Sample_Barcode");
		writer.print('\t');
		writer.print("Matched_Norm_Sample_Barcode");
		writer.print('\t');
		writer.print("Match_Norm_Seq_Allele1");
		writer.print('\t');
		writer.print("Match_Norm_Seq_Allele2");
		writer.print('\t');
		writer.print("Tumor_Validation_Allele1");
		writer.print('\t');
		writer.print("Tumor_Validation_Allele2");
		writer.print('\t');
		writer.print("Match_Norm_Validation_Allele1");
		writer.print('\t');
		writer.print("Match_Norm_Validation_Allele2");
		writer.print('\t');
		writer.print("Verification_Status");
		writer.print('\t');
		writer.print("Validation_Status");
		writer.print('\t');
		writer.print("Mutation_Status");
		writer.print('\t');
		writer.print("Sequencing_Phase");
		writer.print('\t');
		writer.print("Sequence_Source");
		writer.print('\t');
		writer.print("Validation_Method");
		writer.print('\t');
		writer.print("Score");
		writer.print('\t');
		writer.print("BAM_File");
		writer.print('\t');
		writer.print("Sequencer");
		writer.print('\t');
		writer.print("Tumor_Sample_UUID");
		writer.print('\t');
		writer.print("Matched_Norm_Sample_UUID");
		writer.println();
	}
}
