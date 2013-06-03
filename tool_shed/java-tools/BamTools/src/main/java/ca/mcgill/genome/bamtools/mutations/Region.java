package ca.mcgill.genome.bamtools.mutations;

import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.tree.IntervalForest;

public class Region {
	private final String name;
	private final IntervalForest forest;
	private long regionSize = 0;

	public Region(String name) {
		super();
		this.name = name;
		forest = new IntervalForest();
	}
	
	public void add(Marker interval) {
//		Markers markers = forest.query(interval);
//		if(markers.size() > 0) {
//			throw new RuntimeException("Interval already in region: " +interval.getId());
//		}
		regionSize += interval.size();
		forest.add(interval);
	}

	public String getName() {
		return name;
	}

	public long getRegionSize() {
		return regionSize;
	}
	
	public boolean intersects(Marker interval) {
		Markers markers = forest.query(interval);
		return markers.size() > 0;
	}
}
