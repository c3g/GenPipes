package ca.mcgill.genome;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.nio.charset.Charset;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ca.mcgill.genome.parser.FASTQFileParser;
import ca.mcgill.genome.parser.FASTQPEFileParser;
import ca.mcgill.genome.parser.Sequence;
import ca.mcgill.genome.parser.SequenceFileParser;

public class FilterDuplicates {
  private static final Logger log = LoggerFactory.getLogger(FilterDuplicates.class);

  /**
   * @param args
   */
  public static void main(String[] args) {

    DuplicateKmerBuilder duplicateKmerBuilder = new DuplicateKmerBuilder();
    byte qualityOffset = 33;

    Options options = new Options();
    options.addOption("h", "help", false, "print this message");

    Option opt = OptionBuilder.create('i');
    opt.setArgName("input");
    opt.setDescription("Input File");
    opt.setArgs(1);
    opt.setRequired(true);
    options.addOption(opt);

    opt = OptionBuilder.create("i2");
    opt.setArgName("PE input");
    opt.setDescription("PE Input File");
    opt.setArgs(1);
    opt.setRequired(false);
    options.addOption(opt);
    
    options.addOption("p", "prefix", true, "output prefix (default original filename .dup.gz)");
    options.addOption("r", "readOutput", false, "Only output duplicate read names");
    options.addOption("k", "kmer", true, "kmer size (default: "+duplicateKmerBuilder.getKmerSize()+")");
    options.addOption("o", "offset", true, "offset to get kmer (default: "+duplicateKmerBuilder.getOffset()+")");
    options.addOption("Q", "quality", true, "Quality Offset (default: "+qualityOffset+")");

    CommandLineParser parser = new PosixParser();
    try {
      CommandLine line = parser.parse(options, args);
      if (line.hasOption("help")) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp(FilterDuplicates.class.getName(), options);
      } else {
        
        if (line.hasOption('k')) {
          duplicateKmerBuilder.setKmerSize(Integer.parseInt(line.getOptionValue('k')));
        }

        if (line.hasOption('o')) {
          duplicateKmerBuilder.setOffset(Integer.parseInt(line.getOptionValue('o')));
        }

        if (line.hasOption('Q')) {
          qualityOffset = Byte.valueOf(line.getOptionValue('Q'));
        }

        boolean onlyOutputReadNames = line.hasOption('r');
        File inputFile = new File(line.getOptionValue('i'));
        String read1OutputFile;
        if(line.hasOption("p")) {
          read1OutputFile = line.getOptionValue('p')+".read1.gz";
        }
        else {
          read1OutputFile = line.getOptionValue('i')+".dup.read1.gz";
        }
        File inputPEFile = null;
        String read2OutputFile = null;

        SequenceFileParser fileParser;
        long totalFileSize = inputFile.length();
        
        if (line.hasOption("i2")) {
          inputPEFile = new File(line.getOptionValue("i2"));
          if(line.hasOption("p")) {
            read2OutputFile = line.getOptionValue('p')+".read2.gz";
          }
          else {
            read2OutputFile = line.getOptionValue('i')+".dup.read2.gz";
          }

          totalFileSize += inputPEFile.length();
          fileParser = new FASTQPEFileParser(inputFile, inputPEFile, qualityOffset);
          duplicateKmerBuilder.setPairedSequence(true);
        } else {
          fileParser = new FASTQFileParser(inputFile, qualityOffset);
        }

        long lastPercent=-1;
        System.out.println("Build kmer ["+duplicateKmerBuilder.getKmerSize()+"] Hash");
        while (true) {
          Sequence sequence = fileParser.nextSequence();
          if (sequence == null)
            break;
          if(sequence.getSequence().length() <= (duplicateKmerBuilder.getKmerSize() + duplicateKmerBuilder.getOffset())) {
            continue;
          }

          duplicateKmerBuilder.processSequence(sequence);

          long percent = fileParser.getByteCount()*100l/totalFileSize;
          if(percent != lastPercent) {
            lastPercent = percent;
            System.out.print("\rFile read: " + percent + "%");
          }
        }
        System.out.println();
        fileParser.close();

        double nbKmers = duplicateKmerBuilder.getNbKnownSequencesFound().size();
        double total = duplicateKmerBuilder.getTotalNbReads();
        System.out.println("Dup: " + (100.0 - (nbKmers * 100.0d / total)) + '%');        

        // Second pass
        fileParser = new FASTQFileParser(inputFile, qualityOffset);
        totalFileSize = inputFile.length();
        BufferedWriter read1Writer = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(read1OutputFile)), Charset.forName("ASCII")));
        if (inputPEFile != null) {
          FASTQFileParser filePEParser = new FASTQFileParser(inputPEFile, qualityOffset);
          BufferedWriter read2Writer = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(read2OutputFile)), Charset.forName("ASCII")));

          lastPercent=-1;
          System.out.println("Extract 'best' unique reads");
          while (true) {
            Sequence read1Sequence = fileParser.nextSequence();
            Sequence read2Sequence = filePEParser.nextSequence();
            if (read1Sequence == null)
              break;

            String sequence = read1Sequence.getSequence() + read2Sequence.getSequence();
            if(sequence.length() <= (duplicateKmerBuilder.getKmerSize() + duplicateKmerBuilder.getOffset())) {
              if(onlyOutputReadNames) {
                read1Writer.append(read1Sequence.getHeader());
                read1Writer.newLine();
                read2Writer.append(read2Sequence.getHeader());
                read2Writer.newLine();
              }
              continue;
            }

            int qualitySum = 0;
            for(int idx=0; idx < read1Sequence.getBaseQualities().length; idx++) {
              qualitySum += read1Sequence.getBaseQualities()[idx];
            }
            for(int idx=0; idx < read2Sequence.getBaseQualities().length; idx++) {
              qualitySum += read2Sequence.getBaseQualities()[idx];
            }

            boolean isWritable;
            try {
              isWritable = duplicateKmerBuilder.getKmerCountAndMark(sequence, qualitySum / (read1Sequence.getBaseQualities().length + read2Sequence.getBaseQualities().length));
            } catch (SequenceLengthException e) {
              throw new RuntimeException("Kmer too long, or offset to big: " + read1Sequence.getHeader());
            }
            if(isWritable && !onlyOutputReadNames){
              read1Writer.append('@').append(read1Sequence.getHeader());
              read1Writer.newLine();
              read1Writer.append(read1Sequence.getSequence());
              read1Writer.newLine();
              read1Writer.append('+').append(read1Sequence.getHeader());
              read1Writer.newLine();
              for(int idx=0; idx < read1Sequence.getBaseQualities().length; idx++) {
                read1Writer.append((char)(read1Sequence.getBaseQualities()[idx]+33));
              }
              read1Writer.newLine();
              read2Writer.append('@').append(read2Sequence.getHeader());
              read2Writer.newLine();
              read2Writer.append(read2Sequence.getSequence());
              read2Writer.newLine();
              read2Writer.append('+').append(read2Sequence.getHeader());
              read2Writer.newLine();
              for(int idx=0; idx < read2Sequence.getBaseQualities().length; idx++) {
                read2Writer.append((char)(read2Sequence.getBaseQualities()[idx]+33));
              }
              read2Writer.newLine();
            }
            else if(!isWritable && onlyOutputReadNames) {
              read1Writer.append(read1Sequence.getHeader());
              read1Writer.newLine();
              read2Writer.append(read2Sequence.getHeader());
              read2Writer.newLine();
            }

            long percent = fileParser.getByteCount()*100l/totalFileSize;
            if(percent != lastPercent) {
              lastPercent = percent;
              System.out.print("\rFile read: " + percent + "%");
            }
          }
          System.out.println();
          read1Writer.close();
          read2Writer.close();
          fileParser.close();
          filePEParser.close();
        } else {
          lastPercent=-1;
          System.out.println("Extract 'best' unique reads");
          while (true) {
            Sequence read1Sequence = fileParser.nextSequence();
            if (read1Sequence == null)
              break;

            String sequence = read1Sequence.getSequence();
            if(sequence.length() <= (duplicateKmerBuilder.getKmerSize() + duplicateKmerBuilder.getOffset())) {
              if(onlyOutputReadNames) {
                read1Writer.append(read1Sequence.getHeader());
                read1Writer.newLine();
              }
              continue;
            }

            int qualitySum = 0;
            for(int idx=0; idx < read1Sequence.getBaseQualities().length; idx++) {
              qualitySum += read1Sequence.getBaseQualities()[idx];
            }

            boolean isWritable;
            try {
              isWritable = duplicateKmerBuilder.getKmerCountAndMark(sequence, qualitySum / read1Sequence.getBaseQualities().length);
            } catch (SequenceLengthException e) {
              throw new RuntimeException("Kmer too long, or offset to big: " + read1Sequence.getHeader());
            }

            if(isWritable && !onlyOutputReadNames){
              read1Writer.append('@').append(read1Sequence.getHeader());
              read1Writer.newLine();
              read1Writer.append(read1Sequence.getSequence());
              read1Writer.newLine();
              read1Writer.append('+').append(read1Sequence.getHeader());
              read1Writer.newLine();
              for(int idx=0; idx < read1Sequence.getBaseQualities().length; idx++) {
                read1Writer.append((char)(read1Sequence.getBaseQualities()[idx]+33));
              }
              read1Writer.newLine();
            }
            else if(!isWritable && onlyOutputReadNames) {
              read1Writer.append(read1Sequence.getHeader());
              read1Writer.newLine();
            }

            long percent = fileParser.getByteCount()*100l/totalFileSize;
            if(percent != lastPercent) {
              lastPercent = percent;
              System.out.print("\rFile read: " + percent + "%");
            }
          }
          System.out.println();
          read1Writer.close();
          fileParser.close();
        }
      }
    } catch (ParseException exp) {
      log.error("Parsing failed.  Reason: " + exp.getMessage());
    } catch (UnsupportedEncodingException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public static void printGC() {
    runGC();
    System.out.println("Mem: " + (usedMemory() / 1024 / 1024));
  }

  private static void runGC() {
    // It helps to call Runtime.gc()
    // using several method calls:
    for (int r = 0; r < 4; ++r)
      _runGC();
  }

  private static void _runGC() {
    long usedMem1 = usedMemory(), usedMem2 = Long.MAX_VALUE;
    for (int i = 0; (usedMem1 < usedMem2) && (i < 500); ++i) {
      s_runtime.runFinalization();
      s_runtime.gc();
      Thread.yield();

      usedMem2 = usedMem1;
      usedMem1 = usedMemory();
    }
  }

  private static long usedMemory() {
    return s_runtime.totalMemory() - s_runtime.freeMemory();
  }

  private static final Runtime s_runtime = Runtime.getRuntime();
}
