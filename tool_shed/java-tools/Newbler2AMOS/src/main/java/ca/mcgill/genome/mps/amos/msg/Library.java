/**
 * 
 */
package ca.mcgill.genome.mps.amos.msg;

import java.io.PrintStream;

/**
 * 
 */
public class Library implements MessageType {
  private final long iid;
  private final String eid;
  private final Distribution distribution;

  public Library(long iid, String eid, Distribution distribution) {
    this.iid = iid;
    this.eid = eid;
    this.distribution = distribution;
  }

  @Override
  public void write(PrintStream writer) {
    writer.println("{LIB");
    writer.print("iid:");
    writer.println(iid);
    writer.print("eid:");
    writer.println(eid);
    distribution.write(writer);
    writer.println('}');
  }

  public long getIid() {
    return iid;
  }

  public String getEid() {
    return eid;
  }

  public Distribution getDistribution() {
    return distribution;
  }

}
