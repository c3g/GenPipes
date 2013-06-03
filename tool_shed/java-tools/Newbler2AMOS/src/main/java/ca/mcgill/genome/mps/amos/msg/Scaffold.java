/*
 * Copyright 2010© McGill University and Génome Québec Innovation Centre. 
 * All rights reserved.
 *
 */

package ca.mcgill.genome.mps.amos.msg;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 *
 */
public class Scaffold implements MessageType {
  private final long iid;
  private final String eid;
  private final List<Tile> tiles;

  public Scaffold(long iid, String eid) {
    this.iid = iid;
    this.eid = eid;
    this.tiles = new ArrayList<Tile>();
  }

  public List<Tile> getTiles() {
    return tiles;
  }

  public long getIid() {
    return iid;
  }

  public String getEid() {
    return eid;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void write(PrintStream writer) throws IOException {
    writer.println("{SCF");
    writer.print("iid:");
    writer.println(iid);
    writer.print("eid:");
    writer.println(eid);
    for (Tile tile : tiles) {
      tile.write(writer);
    }
    writer.println('}');
  }
}
