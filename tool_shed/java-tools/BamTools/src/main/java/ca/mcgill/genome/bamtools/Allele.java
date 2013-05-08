package ca.mcgill.genome.bamtools;

public class Allele {

	private String order;
	private String strand;
	private Character base;
	public int firstForwardCs = 0;
	public int firstForwardTs = 0;
	public int secondForwardCs = 0;
	public int secondForwardTs = 0;
	public int firstReverseCs = 0;
	public int firstReverseTs = 0;
	public int secondReverseCs = 0;
	public int secondReverseTs = 0;
	
	public Allele() {
		
		// TODO Auto-generated constructor stub
	}
	
	public void setOrder(String order) {
		this.order = order;
	}
	
	public void setStrand(String strand) {
		this.strand = strand;
	}
	
	public String getOrder() {
		return this.order;
	}
	
	public String getStrand() {
		return this.strand;
	}

	public void setBase(Character a) {
		this.base=a;
	}
	public Character getBase() {
		return this.base;
	}

}
