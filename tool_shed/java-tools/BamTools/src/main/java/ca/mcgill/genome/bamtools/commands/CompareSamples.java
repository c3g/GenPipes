package ca.mcgill.genome.bamtools.commands;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ca.mcgill.genome.bamtools.DefaultTool;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.mcb.pcingola.vcf.VcfGenotype;
import ca.mcgill.mcb.pcingola.vcf.VcfHeader;

public class CompareSamples extends DefaultTool {
	private int minDepth = 10;
	private int minGQ = 0;
	private double minQual = 0;
	private boolean testRef = false;
	private boolean onlyHom = false;
	private boolean outputBadCalls = false;
	private String  sampleNames[] = null;

	@Override
	public String getCmdName() {
		return "cmpsamples";
	}

	@Override
	public String getCmdUsage() {
		return super.getCmdUsage();
	}

	private void printUsage(String errMsg) {
		System.out.println(errMsg);
		printUsageHeader();
		System.out.println("\t--vcf             VCFs to load");
//		System.out.println("\t--sampleIndexes   Comma separated list of per vcf index of samples to use. (default: intersect all samples from all vcfs).");
//		System.out.println("\t                  If you submit 3 vcfs and want to compare sample 1 from the first, 3 from the second and 2 from the third input this:");
//		System.out.println("\t                  1,3,2 (no spaces)");
		System.out.println("\t--sampleNames     Comma separated list of per vcf sample names to compare. (default: intersect all samples from all vcfs).");
		System.out.println("\t                  If you submit 3 vcfs and want to compare sample RENAL from the first, KIDNEY from the second and BRAIN from the third input this:");
		System.out.println("\t                  RENAL,KIDNEY,BRAIN (no spaces)");
		System.out.println("\t--mindepth        Minimum depth on which to use a Variant. (default: " + minDepth + ")");
		System.out.println("\t--minQual         Minimum Variant Quality (not genotype quality) on which to use a Variant. (default: " + minQual + ")");
		System.out.println("\t--minGQ           Minimum Genotype Quality (not variant quality) on which to use a Variant. (default: " + minGQ + ")");
		System.out.println("\t--testref         Test the REF at the same position between VCFs. (default: off)");
		System.out.println("\t--onlyhom         Use only Homozygous calls. (default: off)");
		System.out.println("\t--outputBadCalls  Output the calls that fail. (default: off)");
	}

	@Override
	public int run(String[] args) {
		//String sampleIndexStr[] = null;
		try {
			List<File> vcfs = new ArrayList<File>();
	
			for (int idx = 1; idx < args.length; idx++) {
				if (args[idx].equals("--vcf")) {
					idx++;
					vcfs.add(new File(args[idx]));
				} else if (args[idx].equals("--mindepth")) {
					idx++;
					minDepth = Integer.parseInt(args[idx]);
					
				} else if (args[idx].equals("--sampleNames")) {
					idx++;
					sampleNames = Gpr.split(args[idx], ',');
				} else if (args[idx].equals("--testref")) {
					testRef = true;
				} else if (args[idx].equals("--onlyhom")) {
					onlyHom = true;
				} else if (args[idx].equals("--outputBadCalls")) {
					outputBadCalls = true;
				} else if (args[idx].equals("--minQual")) {
					idx++;
					minQual = Integer.parseInt(args[idx]);
				} else if (args[idx].equals("--minGQ")) {
					idx++;
					minGQ = Integer.parseInt(args[idx]);
				}
			}

			if (vcfs.size() == 0) {
				printUsage("vcf not set");
				return 1;
			}
			if(sampleNames != null && sampleNames.length != vcfs.size()) {
				printUsage("If sampleNames is used it must have one value per vcf");
				return 1;
			}
			
//			sampleIndexes = new int[sampleIndexStr.length];
//			for(int idx=0; idx < sampleIndexStr.length; idx++) {
//				sampleIndexes[idx] = Integer.parseInt(sampleIndexStr[idx]);
//				if(sampleIndexes[idx] < 1) {
//					printUsage("sampleIndexes must be >= 1 (one-based)");
//					return 1;
//				}
//			}
			
			compareSnps(vcfs);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}

		return 0;
	}

//	public void compareSnps(File vcf){
//		int nbSnps=0;
//		int nbTests=0;
//		int nbMatch=0;
//
//		VcfFileIterator vcfParser = new VcfFileIterator(vcf.toString());
//		VcfHeader vcfHeader = vcfParser.readHeader();
//		List<String> sampleNames = vcfHeader.getSampleNames();
//		for (VcfEntry vcfEntry : vcfParser) {
//			//Skip indels for now.
//			if(vcfEntry.isInDel() || vcfEntry.getQuality() >= minQual) {
//				continue;
//			}
//
//			nbSnps++;
//			String genotypes[] = new String[sampleNames.size()];
//			boolean isValid=true;
//			for(int idx=0; idx < sampleNames.size(); idx++) {
//				VcfGenotype sampleGenotype = vcfEntry.getVcfGenotype(idx);
//				int depth = Gpr.parseIntSafe(sampleGenotype.get("DP"));
//				if(depth >= minDepth && Gpr.parseIntSafe(sampleGenotype.get("GQ")) >= minGQ) {
//					genotypes[idx] = sampleGenotype.getGenotypeStr();
//				}
//				else {
//					isValid = false;
//					genotypes[idx] = null;
//				}
//			}
//
//			if(isValid) {
//				nbTests++;
//				boolean sameGenotype = true;
//				String genotype = genotypes[0];
//				for(int idx=1; idx < sampleNames.size(); idx++) {
//					if(!genotype.equals(genotypes[idx])) {
//						sameGenotype = false;
//						break;
//					}
//				}
//
//				if(sameGenotype) {
//					nbMatch++;
//				}
//			}
//		}
//		
//		System.out.println("Nb SNPS : " + nbSnps);
//		System.out.println("Nb Tests: " + nbTests);
//		System.out.println("Nb Match: " + nbMatch);
//		System.out.println("% Match : " + (double)nbMatch*100.0/(double)nbTests);
//	}
	

	public void compareSnps(List<File> vcfs){
		Map<String, Map<String,String>> loc2genotype = new HashMap<String, Map<String,String>>();
		Map<String, String> pos2ref = new HashMap<String, String>();
		int totalSamples = 0;
		for(int vcfIdx=0; vcfIdx < vcfs.size(); vcfIdx++) {
			File vcf = vcfs.get(vcfIdx);
			System.err.println("Parsing: "+vcf);
			VcfFileIterator vcfParser = new VcfFileIterator(vcf.toString());
			VcfHeader vcfHeader = vcfParser.readHeader();
			List<String> allSampleNames = vcfHeader.getSampleNames();
			List<Integer> indexesToUse = new ArrayList<Integer>();
			if(sampleNames != null) {
				boolean found = false;
				String sampleName =  sampleNames[vcfIdx];
				for(int idx=0; idx < allSampleNames.size(); idx++) {
					if(sampleName.equals(allSampleNames.get(idx))) {
						indexesToUse.add(idx);
						found = true;
						break;
					}
				}
				if(!found) {
					throw new RuntimeException("Can't find sample: "+sampleName + " in vcf: "+vcf);
				}
			}
			else {
				for(int idx=0; idx < allSampleNames.size(); idx++) {
					indexesToUse.add(idx);
				}
			}
			
			totalSamples += indexesToUse.size();

			for (VcfEntry vcfEntry : vcfParser) {
				//Skip indels for now.
				if(vcfEntry.isInDel() || vcfEntry.getQuality() < minQual) {
					continue;
				}
				

				String key = vcfEntry.getChromosomeName()+':'+vcfEntry.getStart();

				if(testRef) {
					if(!pos2ref.containsKey(key)) {
						pos2ref.put(key, vcfEntry.getRef());
					}
					else {
						if(!vcfEntry.getRef().equals(pos2ref.get(key))) {
							System.err.println("Refs don't match: Original:"+pos2ref.get(key)+" New:"+vcfEntry.getRef()+" NewFile: "+vcf);
						}
					}
				}

				
				for(Integer idx :  indexesToUse) {
					VcfGenotype sampleGenotype;
					try {
						sampleGenotype = vcfEntry.getVcfGenotype(idx);
					}
					catch(RuntimeException e) {
						System.err.println("Problem at: " + vcf + " Pos: " +key);
						throw new RuntimeException(e);
					}
					int depth = Gpr.parseIntSafe(sampleGenotype.get("DP"));
					if(depth >= minDepth && Gpr.parseIntSafe(sampleGenotype.get("GQ")) >= minGQ) {
						if(!onlyHom || sampleGenotype.isHomozygous()) {
							if(!loc2genotype.containsKey(key)) {
								loc2genotype.put(key, new HashMap<String, String>());
							}
	
							loc2genotype.get(key).put(vcf.toString()+'-' + allSampleNames.get(idx), sampleGenotype.getGenotypeStr());
						}
					}
				}
			}
		}
		
		int nbTests=0;
		int nbMatch=0;
		for(String key : loc2genotype.keySet()) {
			Map<String,String> smplGenotype = loc2genotype.get(key);
			if(smplGenotype.size() == totalSamples) {
				nbTests++;
				boolean sameGenotype = true;

				String genotype = null;
				String sampleName = null;
				for(String sample : smplGenotype.keySet()) {
					if(genotype == null) {
						genotype = smplGenotype.get(sample);
						sampleName = sample.substring(0,sample.indexOf('-'));
					}
					else {
						if(!genotype.equals(smplGenotype.get(sample))) {
							sameGenotype = false;
							if(this.outputBadCalls) {
								String testSampleName = sample.substring(0,sample.indexOf('-'));
								System.err.println("Bad: " + key + " Samples: " + sampleName + ',' + testSampleName + " Genotypes:" + genotype + "," + smplGenotype.get(sample));
							}
							break;
						}	
					}
				}

				if(sameGenotype) {
					nbMatch++;
				}
			}
		}

		int nbSnps=loc2genotype.size();
		System.out.println("Nb SNPS : " + nbSnps);
		System.out.println("Nb Tests: " + nbTests);
		System.out.println("Nb Match: " + nbMatch);
		System.out.println("% Match : " + (double)nbMatch*100.0/(double)nbTests);
	}
}
