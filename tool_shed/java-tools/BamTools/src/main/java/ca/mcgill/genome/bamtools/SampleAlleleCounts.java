package ca.mcgill.genome.bamtools;

import gnu.trove.map.hash.TObjectIntHashMap;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.EnumMap;
import java.util.Iterator;
import java.util.concurrent.Callable;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class SampleAlleleCounts implements Closeable {
	public static enum ReadBases {
		A, T, C, G, N
	};
	
	public static enum IndelAtPosition {
		NONE, DEL, INS, BOTH
	};

	public static int DELETION_PREFIX_SIZE=2;

	private final SAMFileReader validationSAMReader;
	private final SAMFileReader referenceSAMReader;
	private final File validationBAMFile;
	private final File referenceBAMFile;
	private final int minMappingQuality;
	private final int minBaseQuality;
	private final IndexedFastaSequenceFile referenceFasta;
	private final ComputeAlleleCounts alleleCountsComputer;

	public SampleAlleleCounts(IndexedFastaSequenceFile referenceFasta, File validationBAMFile, int minMappingQuality, int minBaseQuality) {
		this(referenceFasta, validationBAMFile, null, minMappingQuality, minBaseQuality);
	}

	public SampleAlleleCounts(IndexedFastaSequenceFile referenceFasta, File validationBAMFile, File referenceBAMFile, int minMappingQuality, int minBaseQuality) {
		alleleCountsComputer = new ComputeAlleleCounts();
		this.referenceFasta = referenceFasta;
		this.validationBAMFile = validationBAMFile;
		this.referenceBAMFile = referenceBAMFile;
		this.minMappingQuality = minMappingQuality;
		this.minBaseQuality = minBaseQuality;
		this.validationSAMReader = new SAMFileReader(validationBAMFile);
		validationSAMReader.setValidationStringency(ValidationStringency.SILENT);
		if (!this.validationSAMReader.hasIndex()) { throw new RuntimeException("The bam file " + validationBAMFile + " does not have an index"); }
		if(referenceBAMFile == null) {
			referenceSAMReader = null;
		}
		else {
			referenceSAMReader = new SAMFileReader(referenceBAMFile);
			referenceSAMReader.setValidationStringency(ValidationStringency.SILENT);
			if (!this.referenceSAMReader.hasIndex()) { throw new RuntimeException("The bam file " + referenceBAMFile + " does not have an index"); }
		}
		this.validationSAMReader.enableIndexCaching(true);
		this.validationSAMReader.setValidationStringency(ValidationStringency.SILENT);
	}

	/**
	 * Not thread safe. Run one at a time.
	 * 
	 * @param chromosome
	 * @param position
	 * @return
	 */
	public Callable<AlleleCounts> getCallableAlleleCountsComputer(String chromosome, int position) {
		return getCallableAlleleCountsComputer(chromosome, position, IndelAtPosition.NONE, -1);
	}

	/**
	 * Not thread safe. Run one at a time.
	 * 
	 * @param chromosome
	 * @param position
	 * @return
	 */
	public Callable<AlleleCounts> getCallableAlleleCountsComputer(String chromosome, int position, IndelAtPosition indelAtPos, int indelLength) {
		alleleCountsComputer.setChromosome(chromosome);
		alleleCountsComputer.setPosition(position);
		alleleCountsComputer.setIndelAtPos(indelAtPos);
		alleleCountsComputer.setIndelLength(indelLength);
		
		return alleleCountsComputer;
	}

	@Override
	public void close() throws IOException {
		validationSAMReader.close();
	}

	/**
	 * Not thread safe. Run one at a time.
	 * 
	 * @param chromosome
	 * @param position
	 * @return
	 */
	private TObjectIntHashMap<String> getAlleleCounts(File bamFile, SAMFileReader reader,  String chromosome, int position) {
		TObjectIntHashMap<String> retVal = new TObjectIntHashMap<String>();
		for (ReadBases bases : ReadBases.values()) {
			retVal.put(bases.toString(), 0);
		}

		SAMRecordIterator samIterator = null;
		try {
			samIterator = reader.query(chromosome, position, position, false);
			while (samIterator.hasNext()) {
				SAMRecord record = samIterator.next();
				if (record.getReadUnmappedFlag() || record.getDuplicateReadFlag() || record.getMappingQuality() < minMappingQuality)
					continue;

				for (AlignmentBlock block : record.getAlignmentBlocks()) {
					if (block.getReferenceStart() <= position && position <= (block.getReferenceStart() + block.getLength())) {
						int refOffset = position - block.getReferenceStart();
						int readBaseOffset = block.getReadStart() - 1 + refOffset;
						if (record.getBaseQualities()[readBaseOffset] >= minBaseQuality) {
							byte bases[] = record.getReadBases();
							String baseStr = String.valueOf((char) bases[readBaseOffset]);
							retVal.put(baseStr, (retVal.get(baseStr) + 1));
						}
						break;
					}
				}
			}
		} catch (RuntimeException e) {
			throw new RuntimeException("Problem reading file:" + bamFile, e);
		} finally {
			if(samIterator != null) {
				samIterator.close();
			}
		}
		return retVal;
	}


	/**
	 * Not thread safe. Run one at a time.
	 * 
	 * @param chromosome
	 * @param position
	 * @return
	 */
	private TObjectIntHashMap<String> getIndelAlleleCounts(File bamFile, SAMFileReader reader,  String chromosome, int position, IndelAtPosition indelAtPosition, int indelLength) {
		TObjectIntHashMap<String> retVal = new TObjectIntHashMap<String>();

		SAMRecordIterator samIterator = null;
		try {
			samIterator = reader.query(chromosome, position, position, false);
			EnumMap<CigarOperator, Integer> operatorsToTest = new EnumMap<CigarOperator, Integer>(CigarOperator.class);
			while (samIterator.hasNext()) {
				SAMRecord record = samIterator.next();
				if (record.getReadUnmappedFlag() || record.getDuplicateReadFlag() || record.getMappingQuality() < minMappingQuality)
					continue;
				Cigar cigar = record.getCigar();
				int start = record.getAlignmentBlocks().get(0).getReferenceStart();
				int refOffset = position - start;
				boolean readIndel=false;
				for(CigarElement cigarElement : cigar.getCigarElements()) {
					if(refOffset == -1 && readIndel == false) {
						if(cigarElement.getOperator().equals(CigarOperator.D)) {
							if(!operatorsToTest.containsKey(CigarOperator.D)) {
								operatorsToTest.put(CigarOperator.D, cigarElement.getLength());
							}
							else {
								int length = operatorsToTest.get(CigarOperator.D).intValue();
								if(length != cigarElement.getLength()) {
									System.err.println("Deletion Lengths vary...");
									if(cigarElement.getLength() > length) {
										operatorsToTest.put(CigarOperator.D, cigarElement.getLength());
									}
								}
							}
						}
						else if(cigarElement.getOperator().equals(CigarOperator.I)) {
							if(!operatorsToTest.containsKey(CigarOperator.I)) {
								operatorsToTest.put(CigarOperator.I, cigarElement.getLength());
							}
							else {
								int length = operatorsToTest.get(CigarOperator.I).intValue();
								if(length != cigarElement.getLength()) {
									System.err.println("Insertion Lengths vary...");
									if(cigarElement.getLength() > length) {
										operatorsToTest.put(CigarOperator.I, cigarElement.getLength());
									}
								}
							}
						}
						break;
					}
					else {
						if(cigarElement.getOperator().consumesReferenceBases()) {
							refOffset -= cigarElement.getLength();
						}
					}
				}
			} // while
			samIterator.close();
			if(operatorsToTest.size() > 1) {
				System.err.println("Insertion AND deletion at the same locus: " + chromosome + ':'+position);
			}

			if(indelAtPosition != IndelAtPosition.BOTH) {
				Iterator<CigarOperator> operatorIter = operatorsToTest.keySet().iterator();
				while(operatorIter.hasNext()) {
					CigarOperator operator = operatorIter.next();
					if(indelAtPosition == IndelAtPosition.DEL && operator.equals(CigarOperator.I)) {
						operatorIter.remove();
					}
					else if(indelAtPosition == IndelAtPosition.INS && operator.equals(CigarOperator.D)) {
						operatorIter.remove();
					}
				}
			}

			samIterator = reader.query(chromosome, position, position, false);
			while (samIterator.hasNext()) {
				SAMRecord record = samIterator.next();
				if (record.getReadUnmappedFlag() || record.getDuplicateReadFlag() || record.getMappingQuality() < minMappingQuality)
					continue;
				Cigar cigar = record.getCigar();
				int start = record.getAlignmentBlocks().get(0).getReferenceStart();
				int refOffset = position - start;
				int readOffset=0;
				boolean readIndel=false;
				for(CigarElement cigarElement : cigar.getCigarElements()) {
					if(refOffset == -1 && readIndel == false) {
						if(cigarElement.getOperator().equals(CigarOperator.D)) {
							readIndel = true;
							ReferenceSequence sequence1 = referenceFasta.getSubsequenceAt(chromosome, position, position);
							String bases = new String(sequence1.getBases(), "ASCII");
							if(!retVal.containsKey(bases)) {
								retVal.put(bases, 0);
							}
							retVal.increment(bases);
							break;
						}
						else if(cigarElement.getOperator().equals(CigarOperator.I)) {
							readIndel = true;
							String bases = new String(record.getReadBases(), readOffset-1, cigarElement.getLength()+1, "ASCII");
							if(!retVal.containsKey(bases)) {
								retVal.put(bases, 0);
							}
							retVal.increment(bases);
							break;
						}
						else if(cigarElement.getOperator().equals(CigarOperator.S)) {
							if(!retVal.containsKey("SoftClipped")) {
								retVal.put("SoftClipped", 0);
							}
							retVal.increment("SoftClipped");
						}
						else {
							if(!cigarElement.getOperator().equals(CigarOperator.M)) {
								throw new RuntimeException("Unknown base type: " + cigarElement.getOperator() + " position: " + chromosome + ':' + position);
							}
						}
					}
					else {
						if(cigarElement.getOperator().consumesReferenceBases()) {
							refOffset -= cigarElement.getLength();
						}
						if(cigarElement.getOperator().consumesReadBases()) {
							readOffset += cigarElement.getLength();
						}
					}
				}

				if(!readIndel) {
					if(operatorsToTest.size() > 0) {
						for(CigarOperator operator : operatorsToTest.keySet()) {
							int operatorLength = operatorsToTest.get(operator).intValue();
							if(operator.equals(CigarOperator.D)) {
								readIndel = true;
								ReferenceSequence sequence = referenceFasta.getSubsequenceAt(chromosome, position, position+operatorLength);
								String bases = new String(sequence.getBases(), "ASCII");
								if(!retVal.containsKey(bases)) {
									retVal.put(bases, 0);
								}
								retVal.increment(bases);
							}
							else if(operator.equals(CigarOperator.I)) {
								readIndel = true;
								ReferenceSequence sequence = referenceFasta.getSubsequenceAt(chromosome, position, position);
								String bases = new String(sequence.getBases(), "ASCII");
								if(!retVal.containsKey(bases)) {
									retVal.put(bases, 0);
								}
								retVal.increment(bases);
							}
						}
					}
					else {
						ReferenceSequence sequence = referenceFasta.getSubsequenceAt(chromosome, position, position+indelLength-1);
						String bases = new String(sequence.getBases(), "ASCII");
						if(!retVal.containsKey(bases)) {
							retVal.put(bases, 0);
						}
						retVal.increment(bases);
//						if(!retVal.containsKey("NoIndel")) {
//						retVal.put("NoIndel", 0);
//					}
//					retVal.increment("NoIndel");	
					}
				}
			}
		} catch (RuntimeException e) {
			throw new RuntimeException("Problem reading file:" + bamFile, e);
		} catch (UnsupportedEncodingException e) {
			throw new RuntimeException(e);
		} finally {
			if(samIterator != null) {
				samIterator.close();
			}
		}
		return retVal;
	}

	private class ComputeAlleleCounts implements Callable<AlleleCounts> {
		private String chromosome;
		private int position;
		private IndelAtPosition indelAtPos;
		private int indelLength;

		@Override
		public AlleleCounts call() throws Exception {
			AlleleCounts ac = new AlleleCounts(chromosome, position);
			if(isIndel()) {
				TObjectIntHashMap<String> bamCounts = getIndelAlleleCounts(validationBAMFile, validationSAMReader, chromosome, position, indelAtPos, indelLength);
				TObjectIntHashMap<String> referenceBAMCounts = null;
				if(referenceBAMFile != null) {
					referenceBAMCounts = getIndelAlleleCounts(referenceBAMFile, referenceSAMReader, chromosome, position, indelAtPos, indelLength);
				}
				ac.setBamCounts(bamCounts);
				ac.setReferenceBAMCounts(referenceBAMCounts);
			}
			else {
				TObjectIntHashMap<String> bamCounts = getAlleleCounts(validationBAMFile, validationSAMReader, chromosome, position);
				TObjectIntHashMap<String> referenceBAMCounts = null;
				if(referenceBAMFile != null) {
					referenceBAMCounts = getAlleleCounts(referenceBAMFile, referenceSAMReader, chromosome, position);
				}
				
				ac.setBamCounts(bamCounts);
				ac.setReferenceBAMCounts(referenceBAMCounts);
			}
			return ac;
		}

		public void setChromosome(String chromosome) {
			this.chromosome = chromosome;
		}

		public void setPosition(int position) {
			this.position = position;
		}

		public boolean isIndel() {
			return getIndelAtPos() != IndelAtPosition.NONE;
		}
		public void setIndelLength(int indelLength) {
			this.indelLength = indelLength;
		}

		public IndelAtPosition getIndelAtPos() {
			return indelAtPos;
		}

		public void setIndelAtPos(IndelAtPosition indelAtPos) {
			this.indelAtPos = indelAtPos;
		}
		
	}
}
