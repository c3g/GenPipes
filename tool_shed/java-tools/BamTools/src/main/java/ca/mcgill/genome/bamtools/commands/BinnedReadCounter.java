package ca.mcgill.genome.bamtools.commands;

import gnu.trove.list.array.TDoubleArrayList;

import java.io.File;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.BAMIndexMetaData;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.StringUtil;
import ca.mcgill.genome.bamtools.Bin;
import ca.mcgill.genome.bamtools.DefaultTool;
import ca.mcgill.genome.bamtools.GenomeBins;

public class BinnedReadCounter extends DefaultTool {
	public enum NormalizeBy {CHR, GENOME};
	private int minMapQ = 0;
	private NormalizeBy normalizeBy = NormalizeBy.GENOME;
	private int window = 200;
	private ReferenceSequenceFile reference = null;
	private boolean computeGC = false;

	@Override
	public String getCmdName() {
		return "bincounter";
	}

	@Override
	public String getCmdUsage() {
		return super.getCmdUsage();
	}

	private void printUsage(String errMsg) {
		System.out.println(errMsg);
		printUsageHeader();
		System.out.println("\t--minMapQ   Min mapping quality (default "+minMapQ+")");
		System.out.println("\t--gc        Compute GC% per window (needs --ref)");
		System.out.println("\t--ref       Indexed Reference Genome. (optional)");
		System.out.println("\t--refbam    Reference (or normal in cancer) BAM to count");
		System.out.println("\t--bam       BAM to count");
		System.out.println("\t--norm      Normalize with chromosome: 'chr' or genome: 'genome'. (default "+normalizeBy.toString().toLowerCase()+")");
		System.out.println("\t--window    Window interval (default "+ window +"bp)");
	}

	@Override
	public int run(String[] args) {
		File sampleBAM = null;
		File refBAM = null;

		for (int idx = 1; idx < args.length; idx++) {
			 if (args[idx].equals("--minMapQ")) {
				idx++;
				minMapQ = Integer.parseInt(args[idx]);
			} else if (args[idx].equals("--gc")) {
				computeGC = true;
			} else if (args[idx].equals("--ref")) {
				idx++;
				reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(args[idx]));
			} else if (args[idx].equals("--bam")) {
				idx++;
				sampleBAM = new File(args[idx]);
			} else if (args[idx].equals("--refbam")) {
				idx++;
				refBAM = new File(args[idx]);
			} else if (args[idx].equals("--window")) {
				idx++;
				window = Integer.parseInt(args[idx]);
			} else if (args[idx].equals("--norm")) {
				idx++;
				normalizeBy = Enum.valueOf(NormalizeBy.class, args[idx].toUpperCase());
			}
		}

		if (sampleBAM == null) {
			printUsage("bam not set");
			return 1;
		}
		if(computeGC && reference == null) {
			printUsage("Need a reference to compute GC%");
			return 1;
		}
		
		GenomeBins bamBins = getBins(sampleBAM);
		GenomeBins refBins = getBins(refBAM);
		computeRatio(bamBins, refBins);

		return 0;
	}

	public void computeRatio(GenomeBins bamBins, GenomeBins refBamBins) {
		PrintStream writer = System.out;
		
		writer.print("chr\tstart\tend\tsample_raw\tref_raw\tsample_normalized\tref_normalized\tln(sample/ref)");
		Map<String, TDoubleArrayList> chr2GCBin = new HashMap<String, TDoubleArrayList>();
		if(computeGC) {
			writer.print("\tGC%");
			computeGC(chr2GCBin);
		}
		writer.println();

		for(String chr : bamBins.getChromosomeNames()) {
			TDoubleArrayList gcBins = chr2GCBin.get(chr);
			Bin bins[] = bamBins.getChromosomeBins(chr);
			Bin refBins[] = refBamBins.getChromosomeBins(chr);
			for(int idx=0; idx < bins.length; idx++) {
				double normalizedCount = -1;
				double normalizedRefCount = -1;
				if(normalizeBy == NormalizeBy.GENOME) {
					normalizedCount = (double)bins[idx].getCount()/(double)bamBins.getTotalHits();
					normalizedRefCount = (double)refBins[idx].getCount()/(double)refBamBins.getTotalHits();
				}
				else if(normalizeBy == NormalizeBy.CHR) {
					normalizedCount = (double)bins[idx].getCount()/(double)bamBins.getChromosomeHitCounts(chr);
					normalizedRefCount = (double)refBins[idx].getCount()/(double)refBamBins.getChromosomeHitCounts(chr);
				}
				else {
					throw new RuntimeException("Unknown normalization");
				}

				writer.print(chr);
				writer.print('\t');
				writer.print(bins[idx].getStart());
				writer.print('\t');
				writer.print(bins[idx].getStop());
				writer.print('\t');
				writer.print(bins[idx].getCount());
				writer.print('\t');
				writer.print(refBins[idx].getCount());
				writer.print('\t');
				writer.print(normalizedCount);
				writer.print('\t');
				writer.print(normalizedRefCount);
				writer.print('\t');
				writer.print(Math.log(normalizedCount/normalizedRefCount));
				if(computeGC) {
					writer.print('\t');
					writer.print(gcBins.get(idx));
				}
				writer.println();
			}
		}
	}

	private void computeGC(Map<String, TDoubleArrayList> chr2GCBin) {
		ReferenceSequence ref = null;
        while ((ref = reference.nextSequence()) != null) {
        	final byte[] refBases = ref.getBases();
            StringUtil.toUpperCase(refBases);
        	final int refLength = refBases.length;

        	TDoubleArrayList gcBins = new TDoubleArrayList((int)Math.ceil((double)refLength/(double)window));
        	chr2GCBin.put(ref.getName(), gcBins);

            long totalBases=0;
            long totalGC=0;
            int windowIdx=-1;
            for (int idx=0; idx < refLength; idx++) {
            	if(idx % window == 0) {
            		if(windowIdx != -1) {
            			gcBins.add((double)totalGC/(double)totalBases);
            		}
            		windowIdx++;
            		totalBases=0;
                    totalGC=0;
            	}
            	totalBases++;
            	if(refBases[idx] == 'C' || refBases[idx] == 'G') {
            		totalGC++;
            	}
            }
        	gcBins.add((double)totalGC/(double)totalBases);
        }
	}

	public GenomeBins getBins(File input) {
		System.err.println("Processing: " + input);
		SAMFileReader samReader = new SAMFileReader(input);
		samReader.setValidationStringency(ValidationStringency.SILENT);
		GenomeBins retVal = new GenomeBins(window);

		long totalNbReads = 0;
		for(SAMSequenceRecord seqRecord : samReader.getFileHeader().getSequenceDictionary().getSequences()) {
			retVal.addChromosome(seqRecord.getSequenceName(), seqRecord.getSequenceLength());
			if(samReader.hasIndex()) {
				BAMIndexMetaData metData =  samReader.getIndex().getMetaData(seqRecord.getSequenceIndex());
				totalNbReads += metData.getAlignedRecordCount();
				totalNbReads += metData.getUnalignedRecordCount();
			}
		}
		
		long prevPercent = 0;
		long currentNbReads=0;
		for (SAMRecord record : samReader) {
			if (!record.getReadUnmappedFlag() && record.getMappingQuality() >= minMapQ) {
				retVal.addRecord(record);
			}
			if(totalNbReads != 0) {
				long percent = currentNbReads*100/totalNbReads;
				if(percent != prevPercent) {
					prevPercent = percent;
					System.err.print('\r');
					System.err.print(percent);
				}
				currentNbReads++;
			}
		}
		if(totalNbReads != 0) {
			System.err.println("\r100%");
		}
		samReader.close();
		return retVal;
	}
}
