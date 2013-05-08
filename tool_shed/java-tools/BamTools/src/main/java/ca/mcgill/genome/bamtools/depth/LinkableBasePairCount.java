package ca.mcgill.genome.bamtools.depth;

import gnu.trove.list.TLinkableAdapter;

public class LinkableBasePairCount extends TLinkableAdapter<LinkableBasePairCount> implements BasePairCount {
	private static final long serialVersionUID = 7273875257830849545L;

	private int bpCount;
	private final byte base;
	
	public LinkableBasePairCount(int value, byte base) {
		bpCount = value;
		this.base = base;
	}

	@Override
	public int getBpCount() {
		return bpCount;
	}

	@Override
	public byte getBase() {
		return base;
	}

	@Override
	public void setBpCount(int bpCount) {
		this.bpCount = bpCount;
	}

	@Override
	public void increment() {
		bpCount++;
	}

	@Override
	public void incrementValid() {
		if(bpCount != NO_VALUE) {
			bpCount++;
		}
	}
}
