package ca.mcgill.genome.bamtools.commands;

import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import ca.mcgill.genome.bamtools.DefaultTool;
import ca.mcgill.genome.bamtools.Variant;
import ca.mcgill.genome.bamtools.homology.HomologyComputer;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.mcb.pcingola.vcf.VcfGenotype;
import ca.mcgill.mcb.pcingola.vcf.VcfHeader;

public class BAMHomology extends DefaultTool {
	private byte threads = 1;
	private double hetThreshold = 0.3;
	private int minDepth = 3;
	private boolean printElapsed = false;

	@Override
	public String getCmdName() {
		return "bamhomology";
	}

	@Override
	public String getCmdUsage() {
		return super.getCmdUsage();
	}

	private void printUsage(String errMsg) {
		System.out.println(errMsg);
		printUsageHeader();
		System.out.println("\t--bam           Samples BAM to validate");
		System.out.println("\t--het           At which threshold do we call a valid het. (default: " + hetThreshold + ")");
		System.out.println("\t--mindepth      Minimum depth on which to call a Variant. (default: " + minDepth + ")");
		System.out.println("\t--printElapsed  Prints elapsed percentage");
		System.out.println("\t--threads       Threads to use. (default: "+threads+")");
		System.out.println("\t--vcf           Reference VCF");
	}

	@Override
	public int run(String[] args) {
		File inputVCF = null;
		File bamFile = null;

		for (int idx = 1; idx < args.length; idx++) {
			if (args[idx].equals("--vcf")) {
				idx++;
				inputVCF = new File(args[idx]);
			} else if (args[idx].equals("--threads")) {
				idx++;
				threads = Byte.parseByte(args[idx]);
			} else if (args[idx].equals("--bam")) {
				idx++;
				bamFile = new File(args[idx]);
			} else if (args[idx].equals("--het")) {
				idx++;
				hetThreshold = Double.parseDouble(args[idx]);
			} else if (args[idx].equals("--bam")) {
				idx++;
				bamFile = new File(args[idx]);
			} else if (args[idx].equals("--mindepth")) {
				idx++;
				minDepth = Integer.parseInt(args[idx]);
			} else if (args[idx].equals("--printElapsed")) {
				printElapsed = true;
			}
		}

		if(inputVCF == null) {
			printUsage("vcf not set");
			return 1;
		}
		if(bamFile == null) {
			printUsage("bam not set or is empty");
			return 1;
		}

		computeHomology(inputVCF, bamFile);
		return 0;
	}

	public List<Variant> parseVCF(File vcf) {
		List<Variant> retVal = new ArrayList<Variant>();
		VcfFileIterator vcfParser = new VcfFileIterator(vcf.toString());
		VcfHeader header = vcfParser.readHeader();
		String sample = header.getSampleNames().get(0);
		for (VcfEntry vcfEntry : vcfParser) {
			VcfGenotype vcfGenotype = vcfEntry.getVcfGenotype(0);
			if(!vcfGenotype.isHomozygous()) {
				continue;
			}
			
			Variant variant = new Variant(vcfEntry.getChromosomeName(), vcfEntry.getStart()+1, vcfEntry.getRef(), vcfEntry.getAltsStr());
			int genotype[] = vcfGenotype.getGenotype();
			variant.addSampleGenotype(sample, genotype[0] == 0 ? vcfEntry.getRef() : vcfEntry.getAlts()[genotype[0]-1], genotype[0] == 0 ? vcfEntry.getRef() : vcfEntry.getAlts()[genotype[0]-1], -1, -1);
			retVal.add(variant);
		}
		vcfParser.close();

		return retVal;
	}

	public void computeHomology(File vcf, File bamFile) {
		List<Variant> variantsToTest = parseVCF(vcf);
		Queue<Variant> variants = new ConcurrentLinkedQueue<Variant>(variantsToTest);
		String sampleName = variantsToTest.get(0).getSampleGenotypes().get(0).getSampleName();

		try {
			HomologyComputer homologyComputers[] = new HomologyComputer[threads];
			Thread workers[] = new Thread[threads];
			ThreadGroup homologyGrp = new ThreadGroup("HomologyComputers");
			for(int idx=0; idx < threads; idx++) {
				homologyComputers[idx] = new  HomologyComputer(variants, sampleName, bamFile, hetThreshold, minDepth, 0, 0);
				workers[idx] = new Thread(homologyGrp, homologyComputers[idx]);
				workers[idx].setDaemon(true);
				workers[idx].start();
			}

			long percent = 0;
			int wait = 1000;
			for(Thread worker : workers) {
				while(worker.isAlive()) {
					Thread.sleep(wait);
					if(printElapsed) {
						long currentPercent = (variantsToTest.size()-variants.size())*100l/variantsToTest.size();
						if(percent != currentPercent) {
							percent = currentPercent;
							System.out.print("\rCompletion: " + percent + '%');
						}
					}
				}
			}
			if(printElapsed) {
				System.out.println("\rCompletion: 100%");
			}

			int nbTests = 0;
			int nbPassed = 0;
			for(HomologyComputer homologyComputer : homologyComputers) {
				nbTests += homologyComputer.getNbTests();
				nbPassed += homologyComputer.getNbPassed();
			}
			
			DecimalFormat formatter = new DecimalFormat("#0.00");
			System.out.println("SrcSample,bamFile,NbTests,NbPassed,FractionPassed");
			System.out.print(sampleName);
			System.out.print(',');
			System.out.print(bamFile);
			System.out.print(',');
			System.out.print(nbTests);
			System.out.print(',');
			System.out.print(nbPassed);
			System.out.print(',');
			System.out.println(formatter.format((double)nbPassed*100d/(double)nbTests));
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		} finally {
		}
	}
}
