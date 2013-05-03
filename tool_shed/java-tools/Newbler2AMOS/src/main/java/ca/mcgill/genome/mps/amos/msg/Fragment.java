/**
 * 
 */
package ca.mcgill.genome.mps.amos.msg;

import java.io.PrintStream;
import java.util.EnumMap;
import java.util.Map;

/**
 *
 */
public class Fragment implements MessageType {
  public enum FragmentType {
    Other, BAC, Insert, Transposon, Walk
  };

  private static final Map<FragmentType, Character> fragmentTypeCode;

  private final long iid;
  private final String eid;
  private final Library library;
  private final FragmentType type;
  private final Read readLeft;
  private final Read readRight;

  static {
    fragmentTypeCode = new EnumMap<FragmentType, Character>(FragmentType.class);
    fragmentTypeCode.put(FragmentType.Other, 'X');
    fragmentTypeCode.put(FragmentType.BAC, 'B');
    fragmentTypeCode.put(FragmentType.Insert, 'I');
    fragmentTypeCode.put(FragmentType.Transposon, 'T');
    fragmentTypeCode.put(FragmentType.Walk, 'W');
  };

  public Fragment(long iid, String eid, Library library, FragmentType type, Read read) {
    this(iid, eid, library, type, read, null);
  }

  public Fragment(long iid, String eid, Library library, FragmentType type, Read readLeft, Read readRight) {
    this.iid = iid;
    this.eid = eid;
    this.library = library;
    this.type = type;
    this.readLeft = readLeft;
    this.readRight = readRight;

    if (readLeft != null) {
      readLeft.setFragment(this);
    }
    if (readRight != null) {
      readRight.setFragment(this);
    }
  }

  @Override
  public void write(PrintStream writer) {
    writer.println("{FRG");
    writer.print("iid:");
    writer.println(iid);
    writer.print("eid:");
    writer.println(eid);
    writer.print("lib:");
    writer.println(library.getIid());
    if (readLeft != null && readRight != null) {
      writer.print("rds:");
      writer.print(readLeft.getIid());
      writer.print(',');
      writer.println(readRight.getIid());
    }
    writer.print("typ:");
    writer.println(fragmentTypeCode.get(type));
    writer.println('}');
  }

  public static Map<FragmentType, Character> getFragmenttypecode() {
    return fragmentTypeCode;
  }

  public long getIid() {
    return iid;
  }

  public String getEid() {
    return eid;
  }

  public Library getLibrary() {
    return library;
  }

  public FragmentType getType() {
    return type;
  }

  public Read getRead() {
    return readLeft;
  }

  public Read getReadLeft() {
    return readLeft;
  }

  public Read getReadRight() {
    return readRight;
  }

}
