package ca.mcgill.genome.bamtools;

public class HetSite {

	private final String chromosome;
	private final int position;
	private final char allele1;
	private final char allele2;
	
	public HetSite(String chromosome, int position, char allele1, char allele2) {
		super();
		this.chromosome = chromosome;
		this.position = position;
		this.allele1 = allele1;
		this.allele2 = allele2;
	}
	
	public String getChromosome() {
		return chromosome;
	}
	
	public int getPosition() {
		return position;
	}
	
	public char getAllele1() {
		return allele1;
	}
	public char getAllele2() {
		return allele2;
	}
	
}
