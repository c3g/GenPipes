package ca.mcgill.genome.bamtools.depth;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CoordMath;
import ca.mcgill.genome.bamtools.depth.RegionDepthCounter.RangeFilter;

public class ComputeDepth implements Callable<DepthInterval> {
	private final DepthInterval interval;
	private final File refFasta;
	private final File bam;
	private final boolean ommitN;
	private final boolean computeGC;
	private final int minMappingQuality;
	private final int minBaseQuality;
	
	public ComputeDepth(File bam, File refFasta, DepthInterval interval, boolean ommitN, boolean computeGC, int minMappingQuality, int minBaseQuality) {
		// DepthInterval IS NOT thread-safe. Might need to clone.
		this.bam = bam;
		this.interval = interval;
		this.refFasta = refFasta;
		this.ommitN = ommitN;
		this.computeGC = computeGC;
		this.minMappingQuality = minMappingQuality;
		this.minBaseQuality = minBaseQuality; 
	}

	public void adjustPosition(RegionDepthCounter depthCounter, int position) {
		List<BasePairCount> binList = depthCounter.adjustToPosition(position);
		for(BasePairCount bpCount : binList) {
			if(bpCount.getBpCount() != BasePairCount.NO_VALUE) {
				interval.incrementNbAccessibleBases();
				interval.incrementDepthBin(bpCount.getBpCount());
			}
			
			if(bpCount.getBase() != BasePairCount.MISSING_BASE) {
				interval.incrementNbBases();
				if(bpCount.getBase() == 'G' || bpCount.getBase() == 'C' || bpCount.getBase() == 'g' || bpCount.getBase() == 'c' ) {
					interval.incrementNbGCBases();
				}
			}
		}
	}
	@Override
	public DepthInterval call() throws Exception {
		SAMFileReader reader = null;
		IndexedFastaSequenceFile fastaRef = null;
		
		try {
			if(refFasta != null) {
				fastaRef = new IndexedFastaSequenceFile(refFasta);
			}
			RegionDepthCounter depthCounter = new RegionDepthCounter(fastaRef, interval.getChromosome(), interval.getStart(), ommitN);
			reader = new SAMFileReader(bam);
			reader.setValidationStringency(ValidationStringency.SILENT);
			SAMRecordIterator samIterator = null;
			QualityFilter filter = new QualityFilter();
			try {
				samIterator = reader.query(interval.getChromosome(), interval.getStart(), interval.getEnd(), false);

				// Find the first good record
				SAMRecord record = null;
				while (samIterator.hasNext()) {
					record = samIterator.next();
					if (!record.getReadUnmappedFlag() && !record.getDuplicateReadFlag() && record.getMappingQuality() >= minMappingQuality)
						break;
				}
				if(record == null)
					return interval;
				
				interval.incrementNbReads();
				if(record.getProperPairFlag())
					interval.incrementNbProperlyPaired();

				int start = record.getAlignmentStart();
				if(start < interval.getStart())
					start = interval.getStart();
				byte bases[];
				if(fastaRef != null) {
					bases = depthCounter.bases(interval.getStart(), start);
				}
				else {
					bases = new byte[ start-interval.getStart()+1];
					Arrays.fill(bases, BasePairCount.MISSING_BASE);
				}
				for(int basePos=interval.getStart(); basePos < start; basePos++) {
					byte base = bases[basePos-interval.getStart()];
					if(!ommitN || !depthCounter.ommit(base)) {
						interval.incrementNbAccessibleBases();
						interval.incrementDepthBin(0);
					}

					if(computeGC) {
						interval.incrementNbBases();
						if(base == 'G' || base == 'g' || base == 'C' || base == 'c') {
							interval.incrementNbGCBases();
						}
					}
				}
				samIterator.close();

				samIterator = reader.query(interval.getChromosome(), interval.getStart(), interval.getEnd(), false);
				while (samIterator.hasNext()) {
					record = samIterator.next();
					if (record.getReadUnmappedFlag() || record.getDuplicateReadFlag() || record.getMappingQuality() < minMappingQuality)
						continue;

					interval.incrementNbReads();
					if(record.getProperPairFlag())
						interval.incrementNbProperlyPaired();

					int recStart = record.getAlignmentStart();
					if(recStart < interval.getStart())
						recStart = interval.getStart();
					
					int resize = record.getAlignmentEnd() + record.getAlignmentEnd()-record.getAlignmentStart();
					if(resize > interval.getEnd()) {
						resize = interval.getEnd();
					}
					depthCounter.resizeTo(resize);
					adjustPosition(depthCounter, recStart);

					filter.setQualities(record.getBaseQualities());
					
					for (AlignmentBlock block : record.getAlignmentBlocks()) {
						int blockStart = block.getReferenceStart();
						if(blockStart > interval.getEnd())
							break;

						if(blockStart < interval.getStart())
							blockStart = interval.getStart();

						int blockEnd = CoordMath.getEnd(block.getReferenceStart(), block.getLength());
						if(blockEnd > interval.getEnd())
							blockEnd = interval.getEnd();

						filter.setBlock(block);
						depthCounter.incrementRangeFiltered(blockStart, blockEnd, filter);
					}
				}
				depthCounter.resizeTo(interval.getEnd());
				adjustPosition(depthCounter, interval.getEnd());
			} catch (RuntimeException e) {
				throw new RuntimeException("Problem reading file:" + bam, e);
			} finally {
				if(samIterator != null) {
					samIterator.close();
				}
			}
		}
		finally {
			if(reader != null)
				reader.close();
		}
		return interval;
	}
	
	private class QualityFilter implements RangeFilter {
		private byte qualities[];
		private AlignmentBlock block;

		@Override
		public boolean isValidToIncrement(int rangeIndex, LinkableBasePairCount basePairCount) {
			// getReadStart is one-based like the rest
			return qualities[block.getReadStart()-1+rangeIndex] >= minBaseQuality;
		}

		public void setQualities(byte[] qualities) {
			this.qualities = qualities;
		}

		public void setBlock(AlignmentBlock block) {
			this.block = block;
		}
	}
}
