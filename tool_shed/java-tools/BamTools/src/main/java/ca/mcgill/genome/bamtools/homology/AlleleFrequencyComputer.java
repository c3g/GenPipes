package ca.mcgill.genome.bamtools.homology;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.Queue;
import java.util.concurrent.Callable;

import ca.mcgill.genome.bamtools.AlleleCounts;
import ca.mcgill.genome.bamtools.SampleAlleleCounts;

public class AlleleFrequencyComputer implements Closeable, Runnable {
	private final Queue<AlleleCounts> positions;
	private final SampleAlleleCounts sampleAlleleCounts;
	
	public AlleleFrequencyComputer(Queue<AlleleCounts> positions, File bam, int minMappingQuality, int minBaseQuality) {
		this.positions = positions;
		sampleAlleleCounts = new SampleAlleleCounts(null, bam, minMappingQuality, minBaseQuality);
	}

	@Override
	public void run() {
		for(AlleleCounts alleleCounts = positions.poll(); alleleCounts != null; alleleCounts = positions.poll()) {
			Callable<AlleleCounts> callableAlleleCounts = sampleAlleleCounts.getCallableAlleleCountsComputer(alleleCounts.getChromosome(), (int)alleleCounts.getPosition());
			AlleleCounts alleleCountsValue = null;
			try {
				alleleCountsValue = callableAlleleCounts.call();
			} catch (Exception e) {
				throw new RuntimeException();
			}
			alleleCounts.setBamCounts(alleleCountsValue.getBamCounts());
		}
	}

	@Override
	public void close() throws IOException {
		sampleAlleleCounts.close();
	}
}
