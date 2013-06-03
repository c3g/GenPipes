package ca.mcgill.genome.bamtools;

public class Bin {
	private final String chromosome;
	private final int start;
	private final int stop;
	private int count = 0;

	public Bin(String chromosome, int start, int stop) {
		super();
		this.chromosome = chromosome;
		this.start = start;
		this.stop = stop;
	}

	public void add() {
		count++;
	}

	public String getChromosome() {
		return chromosome;
	}

	public int getStart() {
		return start;
	}

	public int getStop() {
		return stop;
	}

	public int getCount() {
		return count;
	}
	
	
}
