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
public class Contig implements MessageType {
  private final long iid;
  private final String eid;
  private final char status;
  private final String sequence;
  private final byte qualities[];
  private final List<Tile> tiles;

  public Contig(long iid, String eid, char status, String sequence, byte[] qualities) {
    this.iid = iid;
    this.eid = eid;
    this.status = status;
    this.sequence = sequence;
    this.qualities = qualities;
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

  public char getStatus() {
    return status;
  }

  public String getSequence() {
    return sequence;
  }

  public byte[] getQualities() {
    return qualities;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void write(PrintStream writer) throws IOException {
    writer.println("{CTG");
    writer.print("iid:");
    writer.println(iid);
    writer.print("eid:");
    writer.println(eid);
    writer.print("sts:");
    writer.println(status);
    writer.print("seq:");
    int read = 0;
    for (char base : sequence.toCharArray()) {
      if ((read % DEFAULT_LINE_LENGTH) == 0) {
        writer.println();
      }
      read++;
      writer.print(base);
    }
    writer.println();
    writer.println('.');

    writer.print("qlt:");
    read = 0;
    for (byte baseQuality : qualities) {
      if ((read % DEFAULT_LINE_LENGTH) == 0) {
        writer.println();
      }
      read++;
      writer.print((char) (baseQuality + DEFAULT_QUALITY_SCALE));
    }
    writer.println();
    writer.println('.');
    for (Tile tile : tiles) {
      tile.write(writer);
    }
    writer.println('}');
  }
}
