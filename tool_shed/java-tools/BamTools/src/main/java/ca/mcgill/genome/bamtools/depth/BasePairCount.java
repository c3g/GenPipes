package ca.mcgill.genome.bamtools.depth;

public interface BasePairCount {
	public static final int NO_VALUE=-1;
	public static final byte MISSING_BASE=-1;

	public int getBpCount();
	public byte getBase();
	public void setBpCount(int bpCount);
	public void increment();
	public void incrementValid();
}
