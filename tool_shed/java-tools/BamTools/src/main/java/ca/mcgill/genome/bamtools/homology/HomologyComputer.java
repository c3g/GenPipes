package ca.mcgill.genome.bamtools.homology;

import gnu.trove.map.hash.TObjectDoubleHashMap;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.Callable;

import ca.mcgill.genome.bamtools.AlleleCounts;
import ca.mcgill.genome.bamtools.SampleAlleleCounts;
import ca.mcgill.genome.bamtools.Variant;

public class HomologyComputer implements Closeable, Runnable {
	private final Queue<Variant> variants;
	private final SampleAlleleCounts sampleAlleleCounts;
	private final double hetThreshold; 
	private int minDepth;
	private int nbTests;
	private int nbPassed;
	
	public HomologyComputer(Queue<Variant> variants, String sample, File bam, double hetThreshold, int minDepth, int minMappingQuality, int minBaseQuality) {
		this.variants = variants;
		sampleAlleleCounts = new SampleAlleleCounts(null, bam, minMappingQuality, minBaseQuality);
		this.hetThreshold = hetThreshold;
		nbTests = 0;
		nbPassed = 0;
	}

	@Override
	public void run() {
		for(Variant variant = variants.poll(); variant != null; variant = variants.poll()) {
			Callable<AlleleCounts> alleleCountsComputer = sampleAlleleCounts.getCallableAlleleCountsComputer(variant.getChromosome(), variant.getPosition());
			AlleleCounts alleleCounts = null;
			try {
				alleleCounts = alleleCountsComputer.call();
			} catch (Exception e) {
				throw new RuntimeException();
			}
			
			int totalBaseCount = 0;
			for (String base : alleleCounts.getBamCounts().keySet()) {
				totalBaseCount += alleleCounts.getBamCounts().get(base);
			}
			if(totalBaseCount < minDepth) {
				continue;
			}
			
			TObjectDoubleHashMap<String> alleleFrac = AlleleCounts.computeAlleleFractions(alleleCounts.getBamCounts());
			Set<String> keptBases = new HashSet<String>();
			for (String base : alleleFrac.keySet()) {
				if (alleleFrac.get(base) > hetThreshold) {
					keptBases.add(base);
				}
			}
			// There might be more than 3 alleles, with low depth samples, use all of them anyways.
			if(keptBases.size() == 0) {
				System.err.println("No bases to keep: "+variant.getChromosome() + ':' + variant.getPosition());
				continue;
			}
			
			if(variant.getSampleGenotypes().get(0).isHom()) {
				keptBases.remove(variant.getSampleGenotypes().get(0).getAlleleA());
			}
			else {
				keptBases.remove(variant.getSampleGenotypes().get(0).getAlleleA());
				keptBases.remove(variant.getSampleGenotypes().get(0).getAlleleB());
			}

			nbTests++;
			if(keptBases.size() == 0) {
				nbPassed++;
			}
		}
	}

	public int getNbTests() {
		return nbTests;
	}

	public int getNbPassed() {
		return nbPassed;
	}

	@Override
	public void close() throws IOException {
		sampleAlleleCounts.close();
	}
}
