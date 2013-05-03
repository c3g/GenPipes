package ca.mcgill.genome.bamtools.commands;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

import ca.mcgill.genome.bamtools.DefaultTool;
import ca.mcgill.genome.bamtools.maf.MAFWriter;
import ca.mcgill.genome.bamtools.parsers.GeneMapper;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.mcb.pcingola.vcf.VcfGenotype;
import ca.mcgill.mcb.pcingola.vcf.VcfHeader;

public class VCF2MAF extends DefaultTool {
	public static final int NORMAL_VCF_IDX = 0;
	public static final int TUMOR_VCF_IDX = 1;
	
	private File vcfFile;
	private PrintStream output;
	private GeneMapper geneMapper;
	private int minGenotypeQual = 0;
	private int minDepth = 10;
	private double minQual = 10;
	private int minCLR = 45;
	private String tumorName = null;
	
	@Override
	public String getCmdName() {
		return "vcf2maf";
	}

	@Override
	public String getCmdUsage() {
		return super.getCmdUsage();
	}

	private void printUsage(String errMsg) {
		System.out.println(errMsg);
		printUsageHeader();
		System.out.println("\t--vcf             VCFs to load");
		System.out.println("\t--minDepth        Minimum depth on which to use a Variant. (default: " + minDepth + ")");
		System.out.println("\t--minGQ           Minimum tumor genotype quality. (default: " + minGenotypeQual + ")");
		System.out.println("\t--minQual         Minimum variant quality. (default: " + minQual + ")");
		System.out.println("\t--minCLR          Minimum CLR to call somatic. (default: " + minCLR + ")");
		System.out.println("\t--mapper          Genes to keep file. Use to filter for refseq");
		System.out.println("\t--output          Output file");
		System.out.println("\t--forceTumorName  Force Tumor sample name to this value.");
	}
	
	
	@Override
	public int run(String[] args) {
		try {
			for (int idx = 0; idx < args.length; idx++) {
				if (args[idx].equals("--vcf")) {
					idx++;
					vcfFile = new File(args[idx]);
				} else if (args[idx].equals("--output")) {
					idx++;
					if(args[idx].equals("-"))
						output = System.out;
					else
						output = new PrintStream(args[idx], "ASCII");
				} else if (args[idx].equals("--mapper")) {
					idx++;
					geneMapper = new GeneMapper(new File(args[idx]));
				} else if (args[idx].equals("--forceTumorName")) {
					idx++;
					tumorName = args[idx];
				} else if (args[idx].equals("--minGQ")) {
					idx++;
					minGenotypeQual = Integer.parseInt(args[idx]);
				} else if (args[idx].equals("--minQual")) {
					idx++;
					minQual = Double.parseDouble(args[idx]);
				} else if (args[idx].equals("--minDepth")) {
					idx++;
					minDepth = Integer.parseInt(args[idx]);
				} else if (args[idx].equals("--minCLR")) {
					idx++;
					minCLR = Integer.parseInt(args[idx]);
				}
			}
	
			if(vcfFile == null) {
				printUsage("You need to pass --vcf");
				return 1;
			}
			if(output == null) {
				printUsage("You need to pass --output");
				return 1;
			}

			vcf2maf();
		}
		catch(IOException e) {
			throw new RuntimeException(e);
		}
		
		return 0;
	}
	
	public void vcf2maf() {
		VcfFileIterator vcfParser = new VcfFileIterator(vcfFile.toString());
		VcfHeader header = vcfParser.readHeader();
		List<String> sampleNames = header.getSampleNames();
		String normalSample =  sampleNames.get(NORMAL_VCF_IDX);
		String tumorSample =  sampleNames.get(TUMOR_VCF_IDX);
		MAFWriter writer = null;
		try {
			writer = new MAFWriter(output);
			writer.setNormal(normalSample);
			if(tumorName != null)
				writer.setTumor(tumorName);
			else
				writer.setTumor(tumorSample);
			writer.setMapper(geneMapper);
			for (VcfEntry vcfEntry : vcfParser) {
				VcfGenotype tumorGenotype = vcfEntry.getVcfGenotype(TUMOR_VCF_IDX);
				if(	Gpr.parseIntSafe(tumorGenotype.get("GQ")) >= minGenotypeQual
				&&	tumorGenotype.getGenotypeCode() != 0
				&&	vcfEntry.getQuality() >= minQual
				&&	Gpr.parseIntSafe(tumorGenotype.get("DP")) >= minDepth
				&&	vcfEntry.getInfoInt("CLR") >= minCLR)
					writer.write(vcfEntry);
			}
		}	
		finally {
			if(writer != null) {
				try{
					writer.close();
				} catch (IOException e) {
					//Oh Well
					e.printStackTrace();
				}
			}
		}
	}

}
