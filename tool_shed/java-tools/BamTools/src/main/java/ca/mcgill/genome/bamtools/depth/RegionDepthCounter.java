package ca.mcgill.genome.bamtools.depth;

import gnu.trove.list.linked.TLinkedList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;

public class RegionDepthCounter {
	private static final byte n_VALUE = (byte)'n';
	private static final byte N_VALUE = (byte)'N';
	private final TLinkedList<LinkableBasePairCount> regionWindow;
	private final String sequenceName;
	private int regionStartCoordinate;
	private final IndexedFastaSequenceFile fastaRef;
	private final boolean ommitN;
	
	public RegionDepthCounter(IndexedFastaSequenceFile fastaRef, String sequenceName, int startPosition, boolean ommitN) {
		this.fastaRef = fastaRef;
		this.regionWindow = new TLinkedList<LinkableBasePairCount>();
		this.sequenceName=sequenceName;
		this.regionStartCoordinate = startPosition;
		this.ommitN = ommitN;
		
		byte base = BasePairCount.MISSING_BASE;
		if(fastaRef != null) {
			ReferenceSequence sequence = fastaRef.getSubsequenceAt(sequenceName, startPosition, startPosition);
			base = sequence.getBases()[0];
		}

		if(ommitN) {
			if(base == N_VALUE || base == n_VALUE) {
				regionWindow.add(new LinkableBasePairCount(LinkableBasePairCount.NO_VALUE, base));
			}
			else {
				regionWindow.add(new LinkableBasePairCount(0, base));
			}
		}
		else {
			regionWindow.add(new LinkableBasePairCount(0, base));
		}
	}

	public List<BasePairCount> adjustToPosition(int position) {
		List<BasePairCount> retVal = new ArrayList<BasePairCount>();
		int nbToRemove = position - regionStartCoordinate;
		regionStartCoordinate = position;
		if(nbToRemove > 0) {
			for(int idx=0; idx < nbToRemove; idx++) {
				retVal.add(regionWindow.removeFirst());
			}
		}
		return retVal;
	}

	public byte[] bases(int start, int stop) {
		ReferenceSequence sequence = fastaRef.getSubsequenceAt(sequenceName, start, stop);
		byte bases[] = sequence.getBases();
		return bases;
	}

	public boolean ommit(byte base) {
		return base == n_VALUE || base == N_VALUE;
	}

	public void resizeTo(int position) {
		int offset = position - regionStartCoordinate;
		int nbToAdd = offset - regionWindow.size()+1;
		if(nbToAdd > 0) {
			byte bases[];
			if(fastaRef != null) {
				bases = bases(regionStartCoordinate+regionWindow.size()-1, position);
			}
			else {
				bases = new byte[nbToAdd];
				Arrays.fill(bases, BasePairCount.MISSING_BASE);
			}

			for(int idx=0; idx < nbToAdd; idx++) {
				regionWindow.add(new LinkableBasePairCount(ommitN && ommit(bases[idx]) ? LinkableBasePairCount.NO_VALUE : 0, bases[idx]));
			}
		}
	}

	public void incrementPosition(int position) {
		int offset = position - regionStartCoordinate;
		regionWindow.get(offset).incrementValid();
	}

	public void incrementRange(int start, int stop) {
		int offsetStart = start - regionStartCoordinate;
		int offsetStop = stop - regionStartCoordinate;
		
		LinkableBasePairCount tmp = regionWindow.get(offsetStart);
        for (int i = offsetStart; i <= offsetStop; i++) {
        	tmp.incrementValid();
        	tmp = tmp.getNext();
        }
	}
	
	public void incrementRangeFiltered(int start, int stop, RangeFilter filter) {
		int offsetStart = start - regionStartCoordinate;
		int offsetStop = stop - regionStartCoordinate;
		
		LinkableBasePairCount tmp = regionWindow.get(offsetStart);
        for (int i = offsetStart, idx=0; i <= offsetStop; i++,idx++) {
        	if(filter != null && filter.isValidToIncrement(idx, tmp))
        		tmp.incrementValid();
        	tmp = tmp.getNext();
        }
	}

	/**
	 * 
	 * @param position
	 * @return List of passed position values
	 */
	public void ignorePosition(int position) {
		int offset = position - regionStartCoordinate;
		regionWindow.get(offset).setBpCount(LinkableBasePairCount.NO_VALUE);
	}
	
	public static interface RangeFilter {
		public boolean isValidToIncrement(int rangeIndex, LinkableBasePairCount basePairCount);
	}
}
