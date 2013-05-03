/**
 * 
 */
package ca.mcgill.genome.mps.amos.msg;

import java.io.IOException;
import java.io.PrintStream;

/**
 * 
 */
public class Read implements MessageType {
  private final long iid;
  private final String eid;
  private final String sequence;
  private final byte[] qlt;
  private final int clearMin;
  private final int clearMax;
  private Fragment fragment;

  public Read(long iid, String eid, String sequence, byte[] qlt, int clearMin, int clearMax) {
    this.iid = iid;
    this.eid = eid;
    this.sequence = sequence;
    this.qlt = qlt;
    this.clearMin = clearMin;
    this.clearMax = clearMax;
  }

  @Override
  public void write(PrintStream writer) throws IOException {
    if (fragment == null) {
      throw new RuntimeException("Fragment is null on read '" + eid + "'");
    }

    writer.println("{RED");
    writer.print("iid:");
    writer.println(iid);
    writer.print("eid:");
    writer.println(eid);

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
    for (byte baseQuality : qlt) {
      if ((read % DEFAULT_LINE_LENGTH) == 0) {
        writer.println();
      }
      read++;
      writer.print((char) (baseQuality + DEFAULT_QUALITY_SCALE));
    }
    writer.println();
    writer.println('.');
    writer.print("frg:");
    writer.println(fragment.getIid());
    writer.print("clr:");
    writer.print(clearMin);
    writer.print(',');
    writer.println(clearMax);
    writer.println('}');
  }

  public long getIid() {
    return iid;
  }

  public String getEid() {
    return eid;
  }

  public String getSequence() {
    return sequence;
  }

  public byte[] getQlt() {
    return qlt;
  }

  public int getClearMin() {
    return clearMin;
  }

  public int getClearMax() {
    return clearMax;
  }

  public Fragment getFragment() {
    return fragment;
  }

  public void setFragment(Fragment fragment) {
    this.fragment = fragment;
  }

}
