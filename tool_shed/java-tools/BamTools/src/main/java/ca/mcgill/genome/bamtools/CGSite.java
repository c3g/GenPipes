package ca.mcgill.genome.bamtools;

import java.util.HashMap;

public class CGSite {

	private final String chromosome;
	private final int position;
	public Allele allele1;
	public Allele allele2;
	private HashMap<String,Integer> readHashFirstForward;
	private HashMap<String,Integer> readHashFirstReverse;
	private HashMap<String,Integer> readHashSecondForward;
	private HashMap<String,Integer> readHashSecondReverse;
	
	public CGSite(String chromosome, int position) {
		super();
		this.chromosome = chromosome;
		this.position = position;
		this.allele1 = new Allele();
		this.allele2 = new Allele();
		this.readHashFirstForward = new HashMap<String,Integer>();
		this.readHashFirstReverse = new HashMap<String,Integer>();
		this.readHashSecondForward = new HashMap<String,Integer>();
		this.readHashSecondReverse = new HashMap<String,Integer>();
	}
	
	public String getChromosome() {
		return chromosome;
	}
	
	public int getPosition() {
		return position;
	}
	
	public boolean checkRead(String readName, String order, String strand) {
		
		
		if(order.equals("first")) {
			if(strand.equals("forward")) {
				if(this.readHashFirstForward.containsKey(readName))	return true; else return false;
			}
			else if(strand.equals("reverse")) {
				if(this.readHashFirstReverse.containsKey(readName))	return true; else return false;
			}
		}
		else if(order.equals("second")) {
				if(strand.equals("forward")) {
					if(this.readHashSecondForward.containsKey(readName)) return true; else return false;
				}
				else if(strand.equals("reverse")) {
					//System.out.println("ici");
					if(this.readHashSecondReverse.containsKey(readName)) {
						//System.out.println("ici");
						return true;
					}
					else return false;
					}
				}
		return false;
	}
	
	public void addRead(String readName, String order, String strand) {
		if(order.equals("first")) {
			if(strand.equals("forward")) {
				this.readHashFirstForward.put(readName, 1);
			}
			else if(strand.equals("reverse")) {
				this.readHashFirstReverse.put(readName, 1);
			}
		}
		else if(order.equals("second")) {
			if(strand.equals("forward")) {
				this.readHashSecondForward.put(readName,1);
			}
			else if(strand.equals("reverse")) {
				this.readHashSecondReverse.put(readName,1);
			}
		}
	}
	
}
