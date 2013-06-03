/*
 * Copyright 2010© McGill University and Génome Québec Innovation Centre. 
 * All rights reserved.
 *
 */

package ca.mcgill.genome.mps.amos.msg;

import java.io.IOException;
import java.io.PrintStream;

/**
 *
 */
public class ContigEdge implements MessageType {
  public enum NodeAdjacency {
    NORMAL, ANTI_NORMAL, INNIE, OUTIE
  };

  private Contig contig1;
  private Contig contig2;
  private int size;
  private int stdDev;
  private NodeAdjacency adj;
  private final char type = 'M';

  /**
   * {@inheritDoc}
   */
  @Override
  public void write(PrintStream writer) throws IOException {

    writer.println("{CTE");
    writer.print("nds:");
    writer.print(contig1.getIid());
    writer.print(',');
    writer.println(contig2.getIid());
    writer.print("adj:");
    writer.println(adj.toString().charAt(0));
    writer.print("sze:");
    writer.println(size);
    writer.print("std:");
    writer.println(stdDev);
    writer.print("typ:");
    writer.println(type);
    writer.println('}');
  }

}
