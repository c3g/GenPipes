package ca.mcgill.genome.bamtools;

import gnu.trove.map.TObjectLongMap;
import gnu.trove.map.hash.TObjectLongHashMap;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CoordMath;

public class GenomeBins {
	private final List<String> chromosomeNames = new ArrayList<String>();
	private final Map<String, Bin[]> chr2Bin = new HashMap<String, Bin[]>();
	private final TObjectLongMap<String> chrHitCounts = new TObjectLongHashMap<String>();
	private final int window;
	private long totalHits = 0;

	public GenomeBins(int window) {
		this.window = window;

	}

	public void addChromosome(String name, int size) {
		int nbBins = (int) Math.ceil((double) size / (double) window);
		Bin bins[] = new Bin[nbBins];
		chrHitCounts.put(name, 0);
		for (int idx = 0; idx < nbBins; idx++) {
			bins[idx] = new Bin(name, idx * window, idx * window + window - 1);
		}

		chr2Bin.put(name, bins);
		chromosomeNames.add(name);
	}

	public void addRecord(SAMRecord record) {
		String chr = record.getReferenceName();
		for(AlignmentBlock block : record.getAlignmentBlocks()) {
			int start = block.getReferenceStart();
			int end = CoordMath.getEnd(block.getReferenceStart(), block.getLength());
			
			int prevIdx = -1;
			for(int idx = start/window; idx <= (end/window); idx++) {
				if(prevIdx != idx) {
					prevIdx = idx;
					chr2Bin.get(chr)[idx].add();
					chrHitCounts.increment(chr);
					totalHits++;
				}
			}
		}
	}

	public List<String> getChromosomeNames() {
		return chromosomeNames;
	}

	public Bin[] getChromosomeBins(String chromosomeName) {
		return chr2Bin.get(chromosomeName);
	}

	public long getChromosomeHitCounts(String chromosomeName) {
		return chrHitCounts.get(chromosomeName);
	}

	public int getWindow() {
		return window;
	}

	public long getTotalHits() {
		return totalHits;
	}

}
