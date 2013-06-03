package ca.mcgill.genome;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.nio.BufferUnderflowException;
import java.util.EnumMap;
import java.util.concurrent.Callable;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class SampleAlleleCounts implements Closeable {
	public static enum ReadBases {
		A, T, C, G, N
	};

	private final SAMFileReader samReader;
	private final File bamFile;
	private final int minMappingQuality;
	private final int minBaseQuality;
	private final ComputeAlleleCounts alleleCountsComputer;

	public SampleAlleleCounts(String sampleName, File bamFile, int minMappingQuality, int minBaseQuality) {
		alleleCountsComputer = new ComputeAlleleCounts();
		this.bamFile = bamFile;
		this.minMappingQuality = minMappingQuality;
		this.minBaseQuality = minBaseQuality;
		this.samReader = new SAMFileReader(bamFile);
		if (!this.samReader.hasIndex()) { throw new RuntimeException("The bam file " + bamFile + " does not have an index"); }
		this.samReader.enableIndexCaching(true);
		this.samReader.setValidationStringency(ValidationStringency.SILENT);
	}

	/**
	 * Not thread safe. Run one at a time.
	 * 
	 * @param chromosome
	 * @param position
	 * @return
	 */
	public Callable<EnumMap<ReadBases, Integer>> getCallableAlleleCountsComputer(String chromosome, int position) {
		alleleCountsComputer.setChromosome(chromosome);
		alleleCountsComputer.setPosition(position);
		return alleleCountsComputer;
	}

	@Override
	public void close() throws IOException {
		samReader.close();
	}

	/**
	 * Not thread safe. Run one at a time.
	 * 
	 * @param chromosome
	 * @param position
	 * @return
	 */
	public EnumMap<ReadBases, Integer> getAlleleCounts(String chromosome, int position) {
		EnumMap<ReadBases, Integer> retVal = new EnumMap<ReadBases, Integer>(ReadBases.class);
		for (ReadBases bases : ReadBases.values()) {
			retVal.put(bases, 0);
		}

		SAMRecordIterator samIterator = null;
		try {
			samIterator = samReader.query(chromosome, position, position, false);
			while (samIterator.hasNext()) {
				SAMRecord record = samIterator.next();
				if (record.getReadUnmappedFlag() || record.getDuplicateReadFlag() || record.getMappingQuality() < minMappingQuality) continue;

				for (AlignmentBlock block : record.getAlignmentBlocks()) {
					if (block.getReferenceStart() <= position && position <= (block.getReferenceStart() + block.getLength())) {
						int refOffset = position - block.getReferenceStart();
						int readBaseOffset = block.getReadStart() - 1 + refOffset;
						if (record.getBaseQualities()[readBaseOffset] >= minBaseQuality) {
							byte bases[] = record.getReadBases();
							String baseStr = String.valueOf((char) bases[readBaseOffset]);
							ReadBases base = ReadBases.valueOf(baseStr);
							retVal.put(base, (retVal.get(base) + 1));
						}
						break;
					}
				}
			}
		} catch (RuntimeException e) {
			throw new RuntimeException("Problem reading file:" + bamFile, e);
		} finally {
			samIterator.close();
		}
		return retVal;
	}

	private class ComputeAlleleCounts implements Callable<EnumMap<ReadBases, Integer>> {
		private String chromosome;
		private int position;

		@Override
		public EnumMap<ReadBases, Integer> call() throws Exception {
			return getAlleleCounts(chromosome, position);
		}

		public void setChromosome(String chromosome) {
			this.chromosome = chromosome;
		}

		public void setPosition(int position) {
			this.position = position;
		}
	}
}
