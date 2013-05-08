package ca.mcgill.genome.bamtools.commands;

import java.io.File;
import java.text.DecimalFormat;
import java.util.List;

import ca.mcgill.genome.bamtools.DefaultTool;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.mcb.pcingola.vcf.VcfGenotype;
import ca.mcgill.mcb.pcingola.vcf.VcfHeader;

public class GenderEstimator extends DefaultTool {
	private File vcfFile;
	private String sexChr = "X";
	private double maleHomRatio = 0.8;
	private double femaleHomRatio = 0.5;
	private String sample = null;
	private int sampleIdx = -1;
	private int minGenotypeQual = 0;
	private int minDepth = 0;
	private double minQual = 100;
	private char gender = 'A';

	@Override
	public String getCmdName() {
		return "gender";
	}

	@Override
	public String getCmdUsage() {
		return super.getCmdUsage();
	}

	private void printUsage(String errMsg) {
		System.out.println(errMsg);
		printUsageHeader();
		System.out.println("\t--vcf          VCF");
		System.out.println("\t--sexChr       Sex chromosome. (default: " + sexChr + ")");
		System.out.println("\t--femaleRatio  Female Homozygous ratio. Needs to be lower or equal than this. (default: " + femaleHomRatio + ")");
		System.out.println("\t--maleRatio    Male Homozygous ratio. Needs to be higher or equal than this. (default: " + maleHomRatio + ")");
		System.out.println("\t--minDepth     Minimum depth on which to use a Variant. (default: " + minDepth + ")");
		System.out.println("\t--minGQ        Minimum tumor genotype quality. (default: " + minGenotypeQual + ")");
		System.out.println("\t--minQual      Minimum variant quality. (default: " + minQual + ")");
		System.out.println("\t--sample       Sample to evaluate in multi-sample vcf. Default is always the first");
		System.out.println("\t--sampleIdx    Index of sample to evaluate in multi-sample vcf. Default is always the first");
		System.out.println("\t--gender       Gender should be this. Values F,M (optional)");
	}
	
	@Override
	public int run(String[] args) {
		for (int idx = 0; idx < args.length; idx++) {
			if (args[idx].equals("--vcf")) {
				idx++;
				vcfFile = new File(args[idx]);
			} else if (args[idx].equals("--sexChr")) {
				idx++;
				sexChr = Chromosome.simpleName(args[idx]);
			} else if (args[idx].equals("--femaleRatio")) {
				idx++;
				femaleHomRatio = Double.parseDouble(args[idx]);
			} else if (args[idx].equals("--maleRatio")) {
				idx++;
				maleHomRatio = Double.parseDouble(args[idx]);
			} else if (args[idx].equals("--sample")) {
				idx++;
				sampleIdx = -1;
				sample = args[idx];
			} else if (args[idx].equals("--sampleIdx")) {
				idx++;
				sample = null;
				sampleIdx = Integer.parseInt(args[idx]);
			} else if (args[idx].equals("--minGQ")) {
				idx++;
				minGenotypeQual = Integer.parseInt(args[idx]);
			} else if (args[idx].equals("--minQual")) {
				idx++;
				minQual = Double.parseDouble(args[idx]);
			} else if (args[idx].equals("--minDepth")) {
				idx++;
				minDepth = Integer.parseInt(args[idx]);
			} else if (args[idx].equals("--gender")) {
				idx++;
				gender = args[idx].charAt(0);
			}
		}

		if(vcfFile == null) {
			printUsage("You need to pass --vcf");
			return 1;
		}
		if(gender != 'A' && gender != 'F' && gender != 'M') {
			printUsage("gender can only be F or M");
			return 1;
		}

		// Default first sample
		if(sampleIdx == -1 && sample == null) {
			sampleIdx = 0;
		}

		estimateGender();
		
		return 0;
	}
	
	public void estimateGender() {
		VcfFileIterator vcfParser = null;
		try {
			vcfParser = new VcfFileIterator(vcfFile.toString());
			VcfHeader header = vcfParser.readHeader();
			List<String> sampleNames = header.getSampleNames();
			if(sampleIdx == -1 && sample != null) {
				for(int idx=0; idx < sampleNames.size(); idx++) {
					if(sampleNames.get(idx).equals(sample)){
						sampleIdx = idx;
						break;
					}
				}
			}
		
			long totalVariants = 0;
			long totalTestedVariants = 0;
			long totalHomVariants = 0;
			for (VcfEntry vcfEntry : vcfParser) {
				if(!Chromosome.simpleName(vcfEntry.getChromosomeName()).equals(sexChr))
					continue;

				totalVariants++;
				VcfGenotype genotype = vcfEntry.getVcfGenotype(sampleIdx);
				boolean passedGQ = true;
				if(minGenotypeQual > 0) {
					passedGQ = Gpr.parseIntSafe(genotype.get("GQ")) >= minGenotypeQual;
				}
				
				if(	passedGQ
				&&	vcfEntry.getQuality() >= minQual
				&&	Gpr.parseIntSafe(genotype.get("DP")) >= minDepth) {
					totalTestedVariants++;
					if(genotype.isHomozygous())
						totalHomVariants++;
				}
			}
			
			DecimalFormat formatter = new DecimalFormat("#0.00");
			double homRatio = (double)totalHomVariants/(double)totalTestedVariants;
			System.out.println("Nb variants: " + totalVariants);
			System.out.println("Nb tested variants: " + totalTestedVariants);
			System.out.println("Nb hom. variants:   " + totalHomVariants);
			System.out.println("Hom ratio:          " + formatter.format(homRatio)) ;
			System.out.print  ("Gender Estimate is: ") ;
			char genderEst = 'A';
			if(homRatio <= femaleHomRatio) {
				System.out.println("Female");
				genderEst = 'F';
			}
			else if(homRatio >= maleHomRatio) {
				System.out.println("Male");
				genderEst = 'M';
			}
			else {
				System.out.println("Ambiguous");
			}
			
			if(gender != 'A') {
				System.out.print("Test outcome:       ");
				if(gender == genderEst) {
					System.out.println("PASSED");
				}
				else {
					System.out.println("FAILED");
				}
			}
		}
		finally {
			if(vcfParser != null) {
				vcfParser.close();
			}
		}
	}

}
