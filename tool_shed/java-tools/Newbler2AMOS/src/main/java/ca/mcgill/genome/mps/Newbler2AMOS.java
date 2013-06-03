package ca.mcgill.genome.mps;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicLong;

import org.apache.commons.lang.math.IntRange;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.SymbolList;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import au.com.bytecode.opencsv.CSVReader;
import ca.mcgill.genome.mps.amos.msg.Contig;
import ca.mcgill.genome.mps.amos.msg.Distribution;
import ca.mcgill.genome.mps.amos.msg.Fragment;
import ca.mcgill.genome.mps.amos.msg.Library;
import ca.mcgill.genome.mps.amos.msg.Read;
import ca.mcgill.genome.mps.amos.msg.Scaffold;
import ca.mcgill.genome.mps.amos.msg.Tile;
import ca.mcgill.genome.mps.core.parser.ACEParser;
import ca.mcgill.genome.mps.core.parser.SFFFileParser;
import ca.mcgill.genome.mps.core.parser.SoftClippedSequence;
import ca.mcgill.genome.mps.core.parser.TrimStatus;
import ca.mcgill.genome.util.ProgressIndicator;

public class Newbler2AMOS implements Runnable, ProgressIndicator {
  private static final Logger log = LoggerFactory.getLogger(Newbler2AMOS.class);
  private static final String UNMATED_LIBRARY_NAME = "unmated";
  private final PrintStream outputStream;
  private final File assemblyDirectory;
  private final List<SffFileDetail> sffFileDetails;
  private final AtomicLong totalItems = new AtomicLong(0);
  private final AtomicLong processedItems = new AtomicLong(0);
  private Map<String, List<TrimStatus>> trimmedReads;
  private Map<String, List<Read>> fragmentReads;
  private Map<String, Read> reads;
  private Map<String, Library> libraries;
  private Map<String, Fragment> fragments;
  private Map<String, Contig> contigs;
  private Map<Long, Contig> iid2Contig;
  private Map<String, Scaffold> scaffolds;

  public Newbler2AMOS(PrintStream outputStream, File assemblyDirectory, List<SffFileDetail> sffFileDetails) {
    this.outputStream = outputStream;
    this.assemblyDirectory = assemblyDirectory;
    this.sffFileDetails = sffFileDetails;
  }

  @Override
  public void run() {
    totalItems.set(0);
    processedItems.set(0);

    try {
      trimmedReads = parseAssemblyTrimInfo();
      fragmentReads = new HashMap<String, List<Read>>();
      reads = new HashMap<String, Read>();
      libraries = new HashMap<String, Library>();
      fragments = new HashMap<String, Fragment>();
      contigs = new HashMap<String, Contig>();
      iid2Contig = new HashMap<Long, Contig>();
      scaffolds = new HashMap<String, Scaffold>();

      int libraryIID = 1;
      int readIID = 1;
      int fragmentIID = 1;
      Library library = new Library(libraryIID++, UNMATED_LIBRARY_NAME, new Distribution());
      libraries.put(library.getEid(), library);
      for (SffFileDetail sffFileDetail : sffFileDetails) {
        String sffFileName = sffFileDetail.getSffFile().getName();
        String libName = sffFileName.substring(0, sffFileName.length() - 4);
        library = new Library(libraryIID++, libName, new Distribution(sffFileDetail.getInsertMean(), sffFileDetail.getInsertDev()));
        libraries.put(library.getEid(), library);

        SFFFileParser sffParser = new SFFFileParser(sffFileDetail.getSffFile());

        while (true) {
          SoftClippedSequence sequence = (SoftClippedSequence) sffParser.nextSequence();
          if (sequence == null) {
            break;
          }

          String readName = sequence.getHeader();
          if (!trimmedReads.containsKey(readName)) {
            log.warn("TrimStatus doesn't contain read '{}'", readName);
          } else {
            for (TrimStatus readStatus : trimmedReads.get(readName)) {
              // All the positions are
              byte baseQuals[] = Arrays.copyOfRange(sequence.getUnTrimmedBaseQualities(), readStatus.getTrimpoints().getMinimumInteger(), readStatus
                  .getTrimpoints().getMaximumInteger() + 1);
              String seq = sequence.getUnTrimmedSequence().substring(readStatus.getTrimpoints().getMinimumInteger(),
                  readStatus.getTrimpoints().getMaximumInteger() + 1);

              IntRange clearRange = new IntRange(0, readStatus.getTrimpoints().getMaximumInteger() - readStatus.getTrimpoints().getMinimumInteger() + 1);
              Read read;
              if (readStatus.getAccno().endsWith("_left")) {
                seq = revComp(seq);
                byte tmpBaseQuals[] = new byte[baseQuals.length];
                for (int idx = baseQuals.length - 1; idx >= 0; idx--) {
                  tmpBaseQuals[baseQuals.length - 1 - idx] = baseQuals[idx];
                }
                read = new Read(readIID++, readStatus.getAccno(), seq, tmpBaseQuals, clearRange.getMinimumInteger(), clearRange.getMaximumInteger());
              } else {
                read = new Read(readIID++, readStatus.getAccno(), seq, baseQuals, clearRange.getMinimumInteger(), clearRange.getMaximumInteger());
              }

              if (!fragmentReads.containsKey(readName)) {
                fragmentReads.put(readName, new ArrayList<Read>());
              }
              fragmentReads.get(readName).add(read);
              if (reads.put(readStatus.getAccno(), read) != null) {
                log.error("Read already loaded {}", readStatus.getAccno());
              }
            }

            Fragment fragment;
            if (trimmedReads.get(readName).size() > 1) {
              Library lib = libraries.get(readName.subSequence(0, 9));
              fragment = new Fragment(fragmentIID++, readName, lib, Fragment.FragmentType.Insert, fragmentReads.get(readName).get(0), fragmentReads.get(
                  readName).get(1));
            } else {
              Library lib = libraries.get(UNMATED_LIBRARY_NAME);
              fragment = new Fragment(fragmentIID++, readName, lib, Fragment.FragmentType.Insert, fragmentReads.get(readName).get(0));
            }
            fragments.put(fragment.getEid(), fragment);
          }
        }
      } // for details
      outputStream.println("{UNV");
      outputStream.println("eid:afg");
      outputStream.println("com:");
      outputStream.println("generated by /home/lletourn/amos-dev/bin/toAmos");
      outputStream.println("Fri Sep 24 14:08:23 2010");
      outputStream.println('.');
      outputStream.println('}');

      for (Library libraryToWrite : libraries.values()) {
        libraryToWrite.write(outputStream);
      }
      for (Fragment fragment : fragments.values()) {
        fragment.write(outputStream);
      }
      for (List<Read> readsInFragment : fragmentReads.values()) {
        for (Read read : readsInFragment) {
          read.write(outputStream);
        }
      }
      outputStream.flush();

      parseAssemblyContigs();
      for (Contig contig : contigs.values()) {
        contig.write(outputStream);
      }

      parseAssemblyScaffolds();
      for (Scaffold scaffold : scaffolds.values()) {
        scaffold.write(outputStream);
      }
      /*
       * scaffolds = parseContigEdges(); for (Scaffold scaffold : scaffolds.values()) { scaffold.write(outputStream); }
       */
    } catch (IOException e) {
      log.error("Couldn't parse files", e);
    }

  }

  private String revComp(String seq) {
    try {
      // make a DNA SymbolList
      SymbolList symL = DNATools.createDNA(seq);

      // reverse complement it
      symL = DNATools.reverseComplement(symL);

      // prove that it worked
      return symL.seqString().toUpperCase();
    } catch (BioException ex) {
      throw new RuntimeException(ex);
    }
  }

  private void parseAssemblyContigs() {
    File aceFile = new File(assemblyDirectory, "454Contigs.ace");
    ACEParser aceParser = new ACEParser(aceFile);
    List<ca.mcgill.genome.mps.core.parser.ace.Contig> aceContigs = aceParser.getContigs();
    Collections.sort(aceContigs, new Comparator<ca.mcgill.genome.mps.core.parser.ace.Contig>() {
      @Override
      public int compare(ca.mcgill.genome.mps.core.parser.ace.Contig o1, ca.mcgill.genome.mps.core.parser.ace.Contig o2) {
        return o1.getName().compareTo(o2.getName());
      }
    });

    for (ca.mcgill.genome.mps.core.parser.ace.Contig aceContig : aceContigs) {
      String realIID = aceContig.getName().substring("contig".length());
      int contigIID = Integer.parseInt(realIID);
      Contig contig = new Contig(contigIID, aceContig.getName(), 'U', aceContig.getUnTrimmedSequence(), aceContig.getUnTrimmedBaseQualities());
      contigs.put(contig.getEid(), contig);
      iid2Contig.put(contig.getIid(), contig);

      for (ca.mcgill.genome.mps.core.parser.ace.Read aceRead : aceContig.getReads().values()) {
        Read read = reads.get(aceRead.getRealName());

        int clearMin = read.getClearMin();
        int clearMax = read.getClearMax();
        IntRange clear = aceRead.getClear();
        if (clear != null) {
          clearMin = clear.getMinimumInteger();
          clearMax = clear.getMaximumInteger();
        }

        if (aceRead.isComplemented()) {
          int tmp = clearMin;
          clearMin = clearMax;
          clearMax = tmp;
        }

        Tile tile = new Tile(read.getIid(), clearMin, clearMax, aceRead.getOffset(), aceRead.gaps());
        contig.getTiles().add(tile);
      }
    }
  }

  private void parseAssemblyScaffolds() throws IOException {
    CSVReader reader = new CSVReader(new FileReader(new File(assemblyDirectory, "454Scaffolds.txt")), '\t');
    String[] nextLine;

    while ((nextLine = reader.readNext()) != null) {
      // N == gap with specified size
      if (nextLine[4].equals("N")) {
        continue;
      }
      String scaffoldEid = nextLine[0];
      Scaffold scaffold = scaffolds.get(scaffoldEid);
      if (scaffold == null) {
        String iidString = scaffoldEid.substring("scaffold".length());
        long iid = Long.parseLong(iidString);
        scaffold = new Scaffold(iid, scaffoldEid);
        scaffolds.put(scaffold.getEid(), scaffold);
      }

      // W == WGS contig
      // http://www.ncbi.nlm.nih.gov/projects/genome/assembly/agp/AGP_Specification.shtml
      if (nextLine[4].equals("W")) {
        Contig contig = contigs.get(nextLine[5]);
        int offset = Integer.parseInt(nextLine[1]) - 1;
        Tile tile = new Tile(contig.getIid(), 0, (Integer.parseInt(nextLine[2]) - offset), offset, null);
        scaffold.getTiles().add(tile);
      }
    }
    reader.close();
  }

  private Map<String, List<TrimStatus>> parseAssemblyTrimInfo() throws IOException {
    Map<String, List<TrimStatus>> retVal = new HashMap<String, List<TrimStatus>>();
    CSVReader reader = new CSVReader(new FileReader(new File(assemblyDirectory, "454TrimStatus.txt")), '\t');
    String[] nextLine;

    // Header
    reader.readNext();
    while ((nextLine = reader.readNext()) != null) {
      String realAccno = nextLine[0];
      String accno = realAccno;
      if (accno.endsWith("_left")) {
        accno = realAccno.substring(0, realAccno.length() - "_left".length());
      } else if (accno.endsWith("_right")) {
        accno = realAccno.substring(0, realAccno.length() - "_right".length());
      }

      TrimStatus newRead = new TrimStatus(realAccno);
      if (!retVal.containsKey(accno)) {
        retVal.put(accno, new ArrayList<TrimStatus>());
      }
      retVal.get(accno).add(newRead);
      if (retVal.get(accno).size() > 2) {
        log.warn("There are more than 2 fragments for read {}", accno);
      }

      int separator = nextLine[1].indexOf('-');
      if (separator == -1) {
        String errMsg = "Trim range isn't a range '" + accno + "': " + nextLine[1];
        log.error(errMsg);
        throw new RuntimeException(errMsg);
      }
      // They're all 1 based, convert to 0 based
      int trimStart = Integer.parseInt(nextLine[1].substring(0, separator)) - 1;
      int trimStop = Integer.parseInt(nextLine[1].substring(separator + 1)) - 1;
      newRead.setTrimpoints(new IntRange(trimStart, trimStop));
      newRead.setTrimmedLength(new Integer(nextLine[2]));

      separator = nextLine[3].indexOf('-');
      if (separator == -1) {
        String errMsg = "Original trim range isn't a range '" + accno + "': " + nextLine[3];
        log.error(errMsg);
        throw new RuntimeException(errMsg);
      }
      // They're all 1 based, convert to 0 based
      trimStart = Integer.parseInt(nextLine[3].substring(0, separator)) - 1;
      trimStop = Integer.parseInt(nextLine[3].substring(separator)) - 1;
      newRead.setOrigTrimpoints(new IntRange(trimStart, trimStop));
      newRead.setOrigTrimmedLength(new Integer(nextLine[4]));

      newRead.setRawLength(new Integer(nextLine[5]));
    }
    reader.close();

    return retVal;
  }

  @Override
  public long getTotalItems() {
    return totalItems.get();
  }

  @Override
  public long getProcessedItems() {
    return processedItems.get();
  }

  /**
   * @param args
   * @throws IOException
   * @throws InterruptedException
   */
  public static void main(String[] args) throws IOException, InterruptedException {
    File assemblyDir = null;
    PrintStream output = System.out;
    List<SffFileDetail> sffFileDetails = new ArrayList<SffFileDetail>();
    for (int i = 0; i < args.length; i++) {
      if (args[i].equals("--assembly")) {
        i++;
        assemblyDir = new File(args[i]);
      } else if (args[i].equals("-o")) {
        i++;
        if (!args[i].equals("-")) {
          output = new PrintStream(new BufferedOutputStream(new FileOutputStream(args[i])), false, "ASCII");
        }
      } else if (args[i].equals("-f")) {
        i++;
        SffFileDetail sffFileDetail = new SffFileDetail(new File(args[i]));
        if (args.length > i && !args[i + 1].startsWith("-")) {
          i++;
          sffFileDetail.setInsertMean(new BigDecimal(args[i]));
          i++;
          sffFileDetail.setInsertDev(new BigDecimal(args[i]));
        }
        sffFileDetails.add(sffFileDetail);
      }
    }

    boolean error = false;
    for (SffFileDetail sffFileDetail : sffFileDetails) {
      if (!sffFileDetail.getSffFile().isFile()) {
        System.err.println("File '" + sffFileDetail.getSffFile() + "' isn't a file or doesn't exist");
        error = true;
      }
    }

    if (!assemblyDir.isDirectory()) {
      System.err.println("Directory '" + assemblyDir + "' isn't a directory or doesn't exist");
      error = true;
    }

    if (!error) {
      Newbler2AMOS newbler2AMOS = new Newbler2AMOS(output, assemblyDir, sffFileDetails);

      // ThroughputCalculator calc = new ThroughputCalculator(newbler2AMOS);
      // ThroughputDisplay display = new ThroughputDisplay(calc, System.err);
      // calc.setSamplingWindow(30);
      // display.setItemName("variants");
      //
      // calc.start(0);
      // display.start();
      try {
        // Thread t = new Thread(newbler2AMOS, Newbler2AMOS.class.getName());
        // t.run();
        // t.join();
      } finally {
        // display.stop();
        // calc.stop();
        newbler2AMOS.run();
        output.close();
        output = null;
      }
    }
  }

}
