/*
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
    or obtain one at http://www.gnu.org

    contact:
    address: R.E. Voorrips
            Wageningen University and Research Centre - Plant Breeding
            P.O. Box 16
            6700 AA Wageningen
            The Netherlands
    e-mail: roeland.voorrips@wur.nl
*/

package PedigreeSim;

import JSci.maths.statistics.TDistribution; //download from http://www.java2s.com/Code/Jar/j/Downloadjscicore11jar.htm
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

/**
 *
 * @author Roeland Voorrips
 */
public class Main {

    /**
     * @param args the command line arguments
     */

    public static String appName = "PedigreeSim";

    public static void main(String[] args) throws Exception {
        simulate(args);
        //Tools.recombtable();
    } //main

    public static void simulate(String[] args) throws Exception {
        File parFile = getParameterFile(args);
        if (parFile==null) { return; }
        System.out.println("Parameter file="+parFile.getName());
        String err = "";
        HashMap<String, String> parMap = readParameterFile(parFile);
        boolean test = false;
        int testIter = 1;
        boolean bivalBidir = false;
        if (parMap.containsKey(strTest)) {
            test = toBoolean(parMap.get(strTest));
            if (parMap.containsKey(strTestIter))
                testIter = Integer.parseInt(parMap.get(strTestIter));
            if (parMap.containsKey(strBivalentsBidirectional))
                bivalBidir = toBoolean(parMap.get(strBivalentsBidirectional));
        }
        boolean interf = false;
        int ploidy;
        if (parMap==null || !parMap.containsKey(strPloidy)) {
            err = parFile.getName()+": no valid "+strPloidy + " found";
        } else
        if (parMap.containsKey(strMapfunction)) {
            interf = parMap.get(strMapfunction).toUpperCase().equals(strKosambi);
        }
        else {
            if(!parMap.containsKey(strHaplostruct)) {
                //mapfunction needed but not present
                err = parFile.getName()+": no valid "+strMapfunction + " found";
            }
        }
        if (!err.isEmpty()) { 
            System.out.println(err);
            return;
        }
        //all data for PopulationData are available
        ploidy = Integer.parseInt(parMap.get(strPloidy));
        double paralMv = 0.0;
        if (parMap.containsKey(strParalMultivalents)) {
            paralMv = Double.parseDouble(parMap.get(strParalMultivalents));
        }
        double pairedCentro = 0.0;
        if (parMap.containsKey(strPairedCentromeres)) {
            pairedCentro = Double.parseDouble(parMap.get(strPairedCentromeres));
        }
        boolean allowNoChiasmata = true;
        if (parMap.containsKey(strAllowNoChiasmata)) {
            allowNoChiasmata = toBoolean(parMap.get(strAllowNoChiasmata));
        }
        boolean naturalPairing = true;
        if (parMap.containsKey(strNaturalPairing)) {
            naturalPairing = toBoolean(parMap.get(strNaturalPairing));
        }
        String missing = strDefaultMissing;
        if (parMap.containsKey(strMissing)) {
            missing = parMap.get(strMissing);
        }
        int seed;
        if (parMap.containsKey(strRandomseed)) {
            seed =  Integer.parseInt(parMap.get(strRandomseed));
        }
        else seed = new Random().nextInt();
        if (!parMap.containsKey(strChromfile)) {
            err = parFile.getName()+": no valid "+strChromfile + " found";
        }
        if (!parMap.containsKey(strOutput)) {
            err = parFile.getName()+": no valid "+strOutput + " found";
        } else {
            if (!canCreateFile(parMap.get(strOutput)+".t64")) {
                err = parFile.getName()+": cannot create "+strOutput + ".*";
            }
        }
        if (!err.isEmpty()) {
            System.out.println(err);
            return;
        }
        PopulationData popdata = new PopulationData(ploidy, interf,
                paralMv, pairedCentro, naturalPairing, allowNoChiasmata,
                bivalBidir, seed, missing);
        
        try { 
            readChromosomeFile(parMap.get(strChromfile), popdata); 
            if (popdata.chromCount()==0)
                throw new Exception(strChromfile+": no chromosomes");
            if (!popdata.allowNoChiasmata) {
                double minLength = 0.70; //minimum chromosome length if allowNoChiasmata is true
                for (Chromosome chr: popdata.getChromosome()) {
                    if (chr.getLength()<minLength)
                        throw new Exception("if " + strAllowNoChiasmata +
                                " is false the minimum chromosome length is " +
                                (minLength*100.0) + " cM");
                }
            }
        }    
        catch (Exception ex) {
            System.out.println(ex.getMessage());
            return;
        }

        //read or generate&write pedigree:

        ArrayList<String[]> pedList;

        if (test) {
            popdata.testMode = true;
            if (parMap.containsKey(strPrintGametes)) {
                popdata.testPrintGametes = Integer.parseInt(parMap.get(strPrintGametes));
            }
            if (parMap.containsKey(strPrintEachIter)) {
                popdata.testPrintEachIter = toBoolean(parMap.get(strPrintEachIter));
            }
            if (parMap.containsKey(strPrintMapdistances)) {
                popdata.testPrintMapdistances = toBoolean(parMap.get(strPrintMapdistances));
            }
            if (parMap.containsKey(strPrintPubTables)) {
                popdata.testPrintPubTables = toBoolean(parMap.get(strPrintPubTables));
            }
            //create population of just one individual:
            String[] TestParent = new String[] {"Parent", missing, missing};
            pedList = new ArrayList<String[]>();
            pedList.add(TestParent);
            err = popdata.createPopulation(pedList);
            if (parMap.containsKey(strPopsize)) {
                popdata.testMeioseCount =
                        Integer.parseInt(parMap.get(strPopsize));
            }
            else err = strPopsize + " not found in parameter file";
        }
        else {
            if (parMap.containsKey(strPedfile)) {
                //get individuals from pedigree file
                ArrayList<String[]> pedListOrig = readPedigreeFile(parMap.get(strPedfile));
                if (!pedListOrig.get(0)[0].isEmpty()) {
                    System.out.println(""+pedListOrig.get(0)[0]);
                    return;
                }
                pedListOrig.remove(0);
                try {
                    pedList = popdata.sortPedigree(pedListOrig);
                    err = popdata.createPopulation(pedList);
                    //check if any changes were made:
                    int i = pedList.size()-1;
                    while (i>=0 && pedList.get(i)[0].equals(pedListOrig.get(i)[0])) {
                        i--;
                    }
                    if (i>=0) {
                        //changes due to sorting, write the new pedigree to file:
                        writePedigreeFile(parMap.get(strOutput)+".ped",popdata);
                    }
                } catch (Exception ex) {
                    System.out.println("Error in pedigree: "+ex.getMessage());
                    return;
                }
            }
            else {
                if (!parMap.containsKey(strPoptype) ||
                    !parMap.containsKey(strPopsize) ) {
                    err = "Either a valid "+strPedfile+" or both "+strPoptype+" and "+
                            strPopsize+" must be specified";
                }
                else {
                    pedList = generateStandardPopulation(
                            parMap.get(strPoptype).toUpperCase(),
                            Integer.parseInt(parMap.get(strPopsize)),
                            missing);
                    if (pedList.isEmpty()) {
                        err = "Invalid population type or population size==0";
                    }
                    else {
                        err = popdata.createPopulation(pedList);
                        if (err.isEmpty())
                            writePedigreeFile(parMap.get(strOutput)+".ped",popdata);
                    }
                }
            } //if PEDFILE else
        } // if test else

        /* pedigree is now available;
         * read or simulate the HaploStruct
         * and generate output based on map and founder genotypes
         */

        if (err.isEmpty()) {

            if (test) {
                if (parMap.containsKey(strHaplostruct)) {
                    try {
                        readHaploStructFiles(parMap.get(strHaplostruct),popdata);
                    } catch (Exception ex) {
                        System.out.println(ex.getMessage());
                        return;
                    }
                }
                // else
                if (!parMap.containsKey(strMapfile)) {
                    err = strMapfile + " not found in parameter file";
                }
                else {
                    try {
                        readMapFile(parMap.get(strMapfile),popdata);
                    } catch (Exception ex) {
                        err = "Error reading "+parMap.get(strMapfile)+": "+ex.getMessage();
                    }
                }
                if (err.isEmpty()) {
                    popdata.testIter = testIter;
                    testMeiosis(popdata,parMap.get(strOutput),2.63);
                }
            } //if test
            else { //normal simulation

                //read or simulate&write haplostructs:

                if (parMap.containsKey(strHaplostruct)) {
                    try {
                        readHaploStructFiles(parMap.get(strHaplostruct),popdata);
                    } catch (Exception ex) {
                        System.out.println(ex.getMessage());
                        return;
                    }
                }
                else { //no Haplostruct files, simulate haplostrucs
                    simulateHaploStructs(popdata);
                    writeHaploStructFiles(parMap.get(strOutput),popdata);
                }
                if (!err.isEmpty()) {
                    System.out.println(err);
                    return;
                }

                //do we have to produce genotypes?
                if (parMap.containsKey(strMapfile)) {
                    //produce and output genotypes
                    try {
                        readMapFile(parMap.get(strMapfile),popdata);
                    } catch (Exception ex) {
                        err = "Error reading "+parMap.get(strMapfile)+": "+ex.getMessage();
                    }
                    boolean founderGenotypes = false;
                    if (err.isEmpty()) {
                        if (parMap.containsKey(strFounderfile)) {
                            try {
                                readFounderGenotypesFile(parMap.get(strFounderfile),popdata);
                                founderGenotypes = true;
                            } catch (Exception ex) {
                                err = "Error reading "+parMap.get(strFounderfile)+": "+ex.getMessage();
                            }
                        }
                    }
                    if (err.isEmpty()) {
                        try {
                            int whichAllele;
                            if (founderGenotypes) {
                                writeGenotypesFile(parMap.get(strOutput)+"_genotypes.dat",false,popdata);
                                whichAllele = MAX_ALLELE;
                            }
                            else whichAllele = 0;
                            writeGenotypesFile(parMap.get(strOutput)+"_founderalleles.dat",true,popdata);
                            writeAlleledoseFile(parMap.get(strOutput)+"_alleledose.dat",whichAllele,founderGenotypes,popdata);
                        } catch (Exception ex) {
                            err = "Error writing files: "+ex.getMessage();
                        }
                    }
                }
                else {
                    //output possible error:
                    //Haplostruct files given but no Map: nothing to do
                    if ( parMap.containsKey(strHaplostruct) &&
                         !parMap.containsKey(strMapfile) ) {
                        err = parFile.getName()+": when "+strHaplostruct+" is given also "+
                                strMapfile+" must be specified";
                    }
                }
            } // if test else
        }
        
        if (!err.isEmpty()) {
            System.out.println(err);
            return;
        }

    } //simulate

    private static char commentChar = ';';
    //Standard keys:
    private static String strMissing = "MISSING";
    private static String strRandomseed = "SEED";
    private static String strPloidy = "PLOIDY";
    private static String strMapfunction = "MAPFUNCTION";
    private static String strChromfile = "CHROMFILE";
    private static String strPedfile = "PEDFILE";
    private static String strMapfile = "MAPFILE";
    private static String strFounderfile = "FOUNDERFILE";
    private static String strHaplostruct = "HAPLOSTRUCT";
    private static String strOutput = "OUTPUT";
    private static String strPoptype = "POPTYPE";
    private static String strPopsize = "POPSIZE";
    //Advanced keys:
    private static String strParalMultivalents = "PARALLELMULTIVALENTS";
    private static String strPairedCentromeres = "PAIREDCENTROMERES";
    private static String strAllowNoChiasmata = "ALLOWNOCHIASMATA";
    private static String strNaturalPairing = "NATURALPAIRING";
    //Test mode keys:
    private static String strTest = "TEST";
    private static String strTestIter = "TESTITER";
    private static String strBivalentsBidirectional = "BIVALENTSBIDIRECTIONAL";
    private static String strPrintGametes = "PRINTGAMETES";
    private static String strPrintEachIter = "PRINTEACHITER";
    private static String strPrintMapdistances = "PRINTMAPDISTANCES";
    private static String strPrintPubTables = "PRINTPUBTABLES"; //only in combination with specific chrom and map files   
    // Value strings:
    private static String strDefaultMissing = "NA";
    private static String strHaldane = "HALDANE";
    private static String strKosambi = "KOSAMBI";
    private static String strF1 = "F1"; //F1 population from 2 heterozygous parents
    private static String strF2 = "F2"; //F2 population from 2 heterozygous perents
    private static String strBC = "BC"; //BC population from 2 heterozygous parents
    private static String strS1 = "S1"; //progeny from selfing one heterozygous parent
                                       //(equivalent to F2 from two homozygous parents)
    
    public static HashMap<String, String> readParameterFile(File parFile)  {
        if (parFile==null) return null;
        HashMap<String, String> parMap = new HashMap<String, String>();
        BufferedReader in;
        try {
            in = new BufferedReader(new FileReader(parFile));
        } catch (Exception ex) { return null; }
        do {
            String s;
            try {
                s = in.readLine();
            } catch (Exception ex) { s = null; }
            if (s==null) break;
            String[] keyval = getKeyandValue(s);
            if ( ( keyval[0].equals(strPloidy) &&
                    isValidPloidy(keyval[1]) )
                    ||
                 ( keyval[0].equals(strMapfunction) &&
                    (keyval[1].toUpperCase().equals(strHaldane) ||
                     keyval[1].toUpperCase().equals(strKosambi)) )
                     ||
                 ( (keyval[0].equals(strChromfile) ||
                    keyval[0].equals(strPedfile) ||
                    keyval[0].equals(strMapfile) ||
                    keyval[0].equals(strFounderfile)) &&
                      new File(keyval[1]).canRead() )
                     ||
                 ( keyval[0].equals(strHaplostruct) &&
                    (new File(keyval[1]+".hsa").canRead() &&
                     new File(keyval[1]+".hsb").canRead()) )
                     ||
                 ( keyval[0].equals(strOutput) &&
                    !keyval[1].isEmpty() ) //further check needed when (re)writing the files
                     ||
                 ( keyval[0].equals(strPoptype) &&
                    (keyval[1].equals(strF1) ||
                     keyval[1].equals(strF2) ||
                     keyval[1].equals(strBC)) )
                     ||
                 ( (keyval[0].equals(strPopsize) ||
                    keyval[0].equals(strTestIter))&&
                    isPositiveInteger(keyval[1]))
                     ||
                 ( (keyval[0].equals(strAllowNoChiasmata) ||
                    keyval[0].equals(strNaturalPairing) ||
                    keyval[0].equals(strBivalentsBidirectional) ||
                    keyval[0].equals(strTest) ||
                    keyval[0].equals(strPrintEachIter) ||
                    keyval[0].equals(strPrintMapdistances) ||
                    keyval[0].equals(strPrintPubTables) ) &&
                     isBoolean(keyval[1]) )
                     ||
                 ( (keyval[0].equals(strParalMultivalents) ||
                    keyval[0].equals(strPairedCentromeres)) &&
                    isProbability(keyval[1]) )
                     ||   
                 ( keyval[0].equals(strMissing) &&
                    !keyval[1].isEmpty() )
                     ||
                 ( (keyval[0].equals(strRandomseed) ||
                    keyval[0].equals(strPrintGametes) ) &&
                    isInteger(keyval[1]) ) 
               ) {
                        parMap.put(keyval[0], keyval[1]);
            }
        } while(true);
        return parMap;
    } //readParameterFile

    public static void readChromosomeFile(String fName, PopulationData popdata)
            throws FileNotFoundException, Exception {
        BufferedReader in;
        try {
            in = new BufferedReader(new FileReader(new File(fName)));
        } catch(Exception ex) {
            throw new Exception(fName+": cannot open file");
        }
        String s;
        //Find headerline:
        do {
            s=in.readLine();
            if (s==null) break;
            s = s.trim();
        } while (s.isEmpty() || s.charAt(0) == commentChar);
        String[] words = Tools.readWords(s);
        if (words.length<3 || !words[0].toUpperCase().equals("CHROMOSOME") ||
                !words[1].toUpperCase().equals("LENGTH") ||
                !words[2].toUpperCase().equals("CENTROMERE")) {
            throw new Exception(fName+": header line not found");
        }
        //Read data lines:
        int wordsNeeded = popdata.ploidy==2 ? 3 : 5;
        do {
            s=in.readLine(); //skip headerline
            if (s==null) break;
            s = s.trim();
            if (!s.isEmpty() && s.charAt(0) != commentChar) {
                words = Tools.readWords(s);
                if (words.length<wordsNeeded) {
                    throw new Exception(fName+": insufficient data on line '"+s+"'");
                } 
                int i = wordsNeeded-1;
                while (i>=0 && !words[i].equals(popdata.missing)) i--;
                if (i>=0) {
                    throw new Exception(fName+": missing data are not allowed");
                }
                String chromName = words[0];
                double length = Double.parseDouble(words[1])/100.0; //cM to Morgan;
                double centromerePos = Double.parseDouble(words[2])/100.0; //cM to Morgan;
                if (popdata.ploidy==2) {
                    popdata.addChromosome(new Chromosome(chromName,length,
                            centromerePos,popdata));
                }
                else {
                    double prefPairingProb = Double.parseDouble(words[3]);
                    double fracQuadrivalents = Double.parseDouble(words[4]);
                    popdata.addChromosome(new TetraploidChromosome(chromName,length,
                            centromerePos,prefPairingProb,fracQuadrivalents,
                            popdata));
                }
            }
        } while(true);
    } //readChromosomeFile

    public static void readMapFile(String fName, PopulationData popdata)
            throws FileNotFoundException, IOException, Exception {
        BufferedReader in;
        try {
            in = new BufferedReader(new FileReader(new File(fName)));
        } catch(Exception ex) {
            throw new Exception(fName+": cannot open file");
        }
        String s;
        //Find headerline:
        do {
            s=in.readLine();
            if (s==null) break;
            s = s.trim();
        } while (s.isEmpty() || s.charAt(0) == commentChar);
        String[] words = Tools.readWords(s);
        if (words.length<3 || !words[0].toUpperCase().equals("MARKER") ||
                !words[1].toUpperCase().equals("CHROMOSOME") ||
                !words[2].toUpperCase().equals("POSITION")) {
            throw new Exception(fName+": header line not found");
        }
        //Read data lines:
        String lastName="";
        Chromosome chrom = null;
        String[] alleleNames = new String[] {"0","1"};
        Locus locus;
        do {
            s=in.readLine();
            if (s==null) break;
            s = s.trim();
            if (!s.isEmpty() && s.charAt(0) != commentChar) {
                words = Tools.readWords(s);
                if (words.length<3) {
                    throw new Exception("Too few items on line '"+s+"'");
                }
                int i = 2;
                while (i>=0 && !words[i].equals(popdata.missing)) i--;
                if (i>=0) {
                    throw new Exception(fName+": missing data are not allowed");
                }
                String mrkName = words[0];
                String chromName = words[1];
                double pos = Double.parseDouble(words[2])/100.0; //cM to Morgan
                if (!chromName.equals(lastName)) {
                    int chr = 0;
                    do {
                        chrom = popdata.getChrom(chr);
                        chr++;
                    } while(chr<popdata.chromCount() &&
                            !chromName.equals(chrom.getChromName()));
                    if (!chromName.equals(chrom.getChromName())) {
                        throw new Exception("chromName "+chromName+" not present");
                    }
                    lastName = chromName;
                }
                locus = new Locus(mrkName,pos,alleleNames,popdata);
                chrom.addLocus(locus);
            }
        } while(true);
    } //readMapFile

    public static void writePedigreeFile(String fName, PopulationData popdata) throws IOException {
        PrintWriter out = new PrintWriter(
                new BufferedWriter(new FileWriter(fName)));
        out.println("Name\tParent1\tParent2");
        for (Individual ind: popdata.getIndividual()) {
            out.print(ind.getIndivName());
            if (ind.isFounder()) {
                out.println("\t"+popdata.missing+"\t"+popdata.missing);
            }
            else {
                out.println("\t"+ind.getParents()[0].getIndivName()+
                            "\t"+ind.getParents()[1].getIndivName());
            }
        }
        out.flush();
    } //writePedigreeFile

    /**
     * readPedigreeFile
     * @param fName
     * @return an ArrayList<String[] of which the first item is an error
     * @throws IOException Only if readline doesn't work, else error message
     * message ("" if no error) and all following items have 3 Strings:
     * individual name and names of parent1 and parent 2 (parents should be
     * strMissing is individual is founder) 
     */
    public static ArrayList<String[]> readPedigreeFile(String fName)
             throws IOException {
        ArrayList<String[]> result = new ArrayList<String[]>();
        result.add(new String[] {""}); //first item is error message
        BufferedReader in;
        try {
            in = new BufferedReader(new FileReader(new File(fName)));
        } catch(Exception ex) {
            result.get(0)[0] =fName+": cannot open file";
            return result;
        }
        String s;
        //find the header line:
        do {
            s = in.readLine();
            if (s==null) break;
        } while (s.trim().isEmpty() ||
                 s.trim().charAt(0) == commentChar);
        String[] words = Tools.readWords(s);
        if (words.length<3 ||
            !words[0].toUpperCase().equals("NAME") ||
            !words[1].toUpperCase().equals("PARENT1") ||
            !words[2].toUpperCase().equals("PARENT2") ) {
            result.get(0)[0] = fName+": header line not found";
            return result;
        }
        //read the first 3 words of each pedigree line
        do {
            s = in.readLine();
            if (s==null) break;
            s = s.trim();
            if (!s.isEmpty() && s.charAt(0) != commentChar) {
                words = Tools.readWords(s);
                if (words.length<3 ||
                    words[1].charAt(0)==commentChar ||
                    words[2].charAt(0)==commentChar) {
                    result.get(0)[0] = "Incomplete line in "+fName+": "+s;
                    return result;
                }
                result.add(new String[] {words[0], words[1],words[2]});
            }
        } while (true);
        return result;
    }

    public static final char HOMOLOG_SEPARATOR = '_';

    private static boolean testGenotypesfileCaption(String s, int homolog) {
        if (s==null || s.length()<3) return false;
        if (s.charAt(s.length()-2) != HOMOLOG_SEPARATOR) return false;
        try {
            int p = Integer.parseInt(s.substring(s.length()-1));
            if (p != homolog) return false;
            return true;
        } catch (Exception ex) {return false;}
    }

    public static void readFounderGenotypesFile(String fName, PopulationData popdata)
            throws FileNotFoundException, IOException, Exception {
        BufferedReader in;
        try {
            in = new BufferedReader(new FileReader(new File(fName)));
        } catch(Exception ex) {
            throw new Exception(fName+": cannot open file");
        }
        String s;
        //Find headerline:
        do {
            s=in.readLine();
            if (s==null) break;
            s = s.trim();
        } while (s.isEmpty() || s.charAt(0) == commentChar);
        String[] words = Tools.readWords(s);
        if (words.length==0 || !words[0].toUpperCase().equals("MARKER")) {
            throw new Exception("header line not found");
        }
        //read headerline:
        //founders may appear in different order than in popdata.individuals
        if (words.length != 1 + popdata.founderAlleleCount) {
            throw new Exception("must contain one column for the marker name and ploidy columns for each founder");
        }
        if (!words[0].toUpperCase().equals("MARKER")) {
            throw new Exception("Header doesn't start with 'marker'");
        }
        int nf = popdata.founderAlleleCount / popdata.ploidy; //number of founders
        int[] indnr = new int[nf]; //will contain the individual number for each
                                   //founder in order of header
        for (int i=0; i<nf; i++) {
            String iname = words[1+i*popdata.ploidy];
            if (!testGenotypesfileCaption(iname,1)) {
                throw new Exception(iname+" doesn't end on '"+HOMOLOG_SEPARATOR+"1'");
            }
            iname = iname.substring(0,iname.length()-2);
            Individual ind = popdata.getIndiv(iname);
            if (ind==null || ind.getParents()[0]!=null) {
                throw new Exception(iname+" is not a founder");
            }
            indnr[i] = popdata.getIndividual().indexOf(ind);
            for (int hom=1; hom<popdata.ploidy; hom++) { //skip hom==0: already done
                iname = words[1+i*popdata.ploidy+hom];
                if (!testGenotypesfileCaption(iname,hom+1) ||
                    !iname.substring(0,iname.length()-2).equals(ind.getIndivName()) ) {
                    throw new Exception(iname+" found but "+
                            ind.getIndivName()+HOMOLOG_SEPARATOR+(hom+1)+
                            " expected");
                }
            }
        }
        // now all captions seem ok
        // if no duplicates occur all founders are represented:
        for (int i=1; i<nf; i++) {
            for (int j=0; j<i; j++) {
                if (indnr[j]==indnr[i]) {
                    throw new Exception("Founder "+popdata.getIndiv(indnr[i])+
                            " occurs more than once in header");
                }
            }
        }

        //read the lines per marker, in order of map:
        int cix = 0;
        Chromosome chrom = popdata.getChrom(cix);
        int loc=0;
        int wordsNeeded = 1 + popdata.founderAlleleCount;
        do {
            s = in.readLine();
            if (s==null) break;
            s = s.trim();
            if (!s.isEmpty() && s.charAt(0) != commentChar) {
                words = Tools.readWords(s);
                if (words.length != wordsNeeded) {
                   throw new Exception("not all lines have "+
                           wordsNeeded+" items");
                }
                int i = wordsNeeded -1;
                while (i>=0 && !words[i].equals(popdata.missing)) i--;
                if (i>=0) {
                    throw new Exception("missing data are not allowed");
                }
                String mrkName = words[0];
                if (chrom==null ||
                    loc>=chrom.getLocus().size() ||
                    !chrom.getLocus().get(loc).getLocusName().equals(mrkName)) {
                    throw new Exception("invalid marker: "+words[0]);
                }
                String[] alleleNames = new String[popdata.founderAlleleCount];
                for (i=0; i<nf; i++) {
                    int startAllele = popdata.getIndiv(indnr[i]).
                            getHaploStruct(cix,0).
                            getFounderAt(chrom.getStartPos());
                    for (int hom=0; hom<popdata.ploidy; hom++) {
                        alleleNames[startAllele+hom] = words[1+i*popdata.ploidy+hom];
                    }
                }
                chrom.getLocus().get(loc).setAlleleName(alleleNames);
                loc++;
                if (loc==chrom.getLocus().size()) {
                    loc = 0;
                    cix++;
                    if (cix==popdata.chromCount()) chrom=null;
                    else chrom = popdata.getChrom(cix);
                }
            }
        } while(true);
        if (chrom!=null) {
            throw new Exception(fName+": end of file found while reading chromosome "+
                    chrom.getChromName());
        }
    } //readFounderGenotypesFile

    /**
     * writeGenotypesFile writes a tab-separated text file with for each 
     * individual and each locus on the map the popdata.ploidy alleles.
     * @param fName the file name; an existing file will be overwritten
     * @param founders if true the founder alleles are written, else the
     * corresponding allele names
     * @param popdata
     * @throws IOException
     * @throws Exception 
     */
    public static void writeGenotypesFile(String fName, boolean founders, PopulationData popdata)
            throws IOException, Exception {
        System.out.println("write file "+fName);
        //File outFile = new File(fName);
        //PrintWriter out = new PrintWriter(new FileWriter(outFile));
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(fName)));
        // write header line:
        out.print("marker");
        for (int i=0; i<popdata.indivCount(); i++) {
            for (int h=0; h<popdata.ploidy; h++) {
                out.print("\t"+ popdata.getIndiv(i).getIndivName()+"_"+(h+1));
            }
        }
        out.println();
        for (int c=0; c<popdata.chromCount(); c++) {
            Chromosome chrom = popdata.getChrom(c);
            System.out.println("chrom "+chrom.getChromName());
            ArrayList<Locus> loci = chrom.getLocus();
            for (int loc=0; loc<loci.size(); loc++) {
                //instead of building a line as a string in memeory and printing
                //it at once it is much faster to print every element separately:
                out.print(loci.get(loc).getLocusName());
                for (int i=0; i<popdata.indivCount(); i++) {
                    if (founders) {
                        for (int h=0; h<popdata.ploidy; h++) {
                            out.print("\t"+ popdata.getIndiv(i).getHaploStruct(c,h).
                                    getFounderAt(loci.get(loc).getPosition()));
                        }
                    } else { //not founders but locus alleles
                        String[] indall = popdata.getIndiv(i).getLocusAllele(c, loc);
                        for (int h=0; h<popdata.ploidy; h++) {
                            out.print("\t"+ indall[h]);
                        }
                    }
                }
                out.println();
            }
        }
        out.flush();
    } //writeGenotypesFile

    public static final int MIN_ALLELE=0;
    public static final int MAX_ALLELE=1;

    /**
     * writeAlleleDoseFile:
     * writes a tab-separated text file with for each individual and each
     * locus on the map the dosage (number of copies) of the allele specified
     * in whichAllele
     * @param fName the file name; an existing file will be overwritten
     * @param whichAllele specifies for which of the alleles the dosage will
     * be written. If actualAlleles is true  and the value if whichAllele
     * is MIN_ALLELE or MAX_ALLELE the allele with
     * the lowest or highest allele name will be counted, else the allele
     * corresponding with founder allele 0 will be counted. If actualAlleles
     * is false the founder allele corresponding to whichAllele is counted
     * (and if whichAllele is outside 0 ... popdata.founderAlleleCount-1
     * the counts will all be 0)
     * @param actualAlleles if true the counting is based on the actual
     * alleles corresponding to each founder allele; if false it is based
     * on the founder alleles directly (and then whichAllele specifies
     * the founder allele counted)
     * @param popdata
     * @throws IOException
     * @throws Exception 
     */
    public static void writeAlleledoseFile(String fName, int whichAllele,
            boolean actualAlleles, PopulationData popdata)
            throws IOException, Exception {
        System.out.println("write file "+fName);
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(fName)));
        // write header line:
        out.print("marker");
        for (int i=0; i<popdata.indivCount(); i++) {
            out.print("\t"+ popdata.getIndiv(i).getIndivName());
        }
        out.println();
        for (int c=0; c<popdata.chromCount(); c++) {
            Chromosome chrom = popdata.getChrom(c);
            System.out.println("chrom "+chrom.getChromName());
            ArrayList<Locus> loci = chrom.getLocus();
            for (int loc=0; loc<loci.size(); loc++) {
                out.print(loci.get(loc).getLocusName());
                double locpos = loci.get(loc).getPosition();
                String matchName="";
                if (actualAlleles) {
                    switch (whichAllele) {
                        case MIN_ALLELE:
                             matchName = loci.get(loc).getMinAlleleName(); break;
                        case MAX_ALLELE:
                             matchName = loci.get(loc).getMaxAlleleName(); break;
                        default: //FOUNDER0_ALLELE
                             matchName = loci.get(loc).getAlleleName(0); break;
                    }
                }
                for (int i=0; i<popdata.indivCount(); i++) {
                    int dose = 0;
                    if (actualAlleles) {
                        String[] indall = popdata.getIndiv(i).getLocusAllele(c, loc);
                        for (int h=0; h<popdata.ploidy; h++) {
                            if (indall[h].equals(matchName)) dose++;
                        }
                    }
                    else { // count founder allele 0 instead of actual allele
                        for (int h=0; h<popdata.ploidy; h++) {
                            if (popdata.getIndiv(i).getHaploStruct(c,h).getFounderAt(locpos)
                                    == whichAllele) dose++;
                        }
                    }
                    out.print("\t" + dose);
                }
                out.println();
            }
        }
        out.flush();
    } //writeAlleledoseFile

    /**
     * 2 files store the HaploStructs of the whole population:
     * *.hsa (the founder segments: set of integers) and
     * *.hsb (the recombination points: set of doubles)
     * for each chromosome the max. length of each haplostruct is used
     * for the whole file
     * @param fName
     * @param popdata
     * @throws IOException
     */
    public static void writeHaploStructFiles(String fName, PopulationData popdata)
            throws IOException {
        PrintWriter hsa = new PrintWriter(
                new BufferedWriter(new FileWriter(fName+".hsa")));
        PrintWriter hsb = new PrintWriter(
                new BufferedWriter(new FileWriter(fName+".hsb")));
        int segcount = 0;
        for (Chromosome chrom: popdata.getChromosome()) {
            int n = chrom.maxSegCount();
            if (n>segcount) segcount = n;
        }

        for (Individual ind: popdata.getIndividual()) {
            for (Chromosome chrom: popdata.getChromosome()) {
                int chrNum = chrom.getChromNumber();
                for (int p=0; p<popdata.ploidy; p++) {
                    hsa.print(ind.getIndivName()+"\t"+chrom.getChromName()+"\t"+(p+1));
                    HaploStruct hs = ind.getHaploStruct(chrNum,p);
                    int hslen = hs.segmentCount();
                    for (int s=0; s<hslen; s++) {
                        hsa.print("\t"+hs.getFounder().get(s));
                        if (s>0) {
                            hsb.print("\t" + hs.getRecombPos().get(s)*100.0); //Morgan to cM
                        }
                    }
                    for (int s=hslen; s<segcount; s++) {
                        hsa.print("\t"+popdata.missing);
                        hsb.print("\t"+popdata.missing);
                    }
                    hsa.println();
                    hsb.println();
                }
            }
        }
        hsa.flush();
        hsb.flush();
    }

    public static void readHaploStructFiles(String fName, PopulationData popdata)
            throws Exception {
        //popdata.clearAllHaploStruct();
        //advance information: chromosomes(name, length etc; ploidy; pedigree
        // hsa has the info on individual and chromosome, and the sequences of founder alleles:
        BufferedReader hsa = new BufferedReader(new FileReader(new File(fName+".hsa")));
        // hsb has the sequences of recombination positions:
        BufferedReader hsb = new BufferedReader(new FileReader(new File(fName+".hsb")));
        //first we make sure that all current haplotypes of individuals are removed
        //to allow a check for individuals occurring multiple times in the haplofiles
        //(Note that founders already get haplotypes in their constructor)
        for (Individual ind: popdata.getIndividual()) {
            ind.clearAllHaplostruct();
        }
        String s;
        int line = 0;
        String iName = "";
        HaploStruct[][] hs = null;
        do {
            s=hsa.readLine();
            if (s!=null) s = s.trim();
            if (s==null || s.equals("")) break;
            String[] words = Tools.readWords(s);
            if (words.length<4) {
                throw new Exception(fName+".hsa: line "+(line+1)+" has too few elements");
            }
            if (line %(popdata.chromCount()*popdata.ploidy) == 0) {
                Individual ind = popdata.getIndiv(words[0]);
                if (ind==null) {
                    throw new Exception(fName+".hsa: individual "+words[0]+" doesn't exist");
                }
                if (ind.getAllHaploStruct()!=null) {
                    throw new Exception(fName+".hsa: individual "+words[0]+" occurs more than once");
                }
                iName = words[0];
                hs = new HaploStruct[popdata.chromCount()][popdata.ploidy];
            }
            else if (!words[0].equals(iName)) {
                throw new Exception(fName+".hsa/hsb: Insufficient lines for individual "+iName);
            }
            int indline = line % (popdata.chromCount()*popdata.ploidy); //nth line for this individual
            int chrnum = indline / popdata.ploidy;                      //chromosome number
            int homol = indline % popdata.ploidy;                       //homolog number
            Chromosome chrom = popdata.getChrom(chrnum);
            if (!words[1].equals(popdata.getChrom(chrnum).getChromName())) {
                throw new Exception(fName+".hsa/hsb, line "+(line+1)+": chromosome "+
                        chrom.getChromName()+" expected but "+
                        words[1]+" found");
            }
            else {
                int homnr;
                try { homnr = Integer.valueOf(words[2]);
                } catch (Exception ex) {
                    throw new Exception(fName+".hsa/hsb, line "+(line+1)+": homolog "+
                        (homol+1) + " expected, but "+words[2]+" found");
                }
                if (homnr != (homol+1)) {
                    throw new Exception(fName+".hsa/hsb, line "+(line+1)+": homolog "+
                            (homol+1) + " expected, but "+words[2]+" found");
                }
            }
            String hsbline = hsb.readLine();
            if (hsbline==null) {
                throw new Exception(fName+".hsb: less lines than "+fName+".hsa");
            }
            int founderAllele = -1;
            try {
                founderAllele = Integer.valueOf(words[3]);
            } catch (Exception ex) {}
            if (founderAllele<0 || founderAllele>=popdata.founderAlleleCount ) {
                throw new Exception(fName+".hsa, line "+(line+1)+": invalid founder allele '"+words[3]+"'");
            }
            String[] hsbwords = Tools.readWords(hsbline);
            if (hsbwords.length != words.length-4) {
                throw new Exception(fName+".hsa/hsb, line "+(line+1)+": number of items in both files don't match");
            }
            HaploStruct hast = new HaploStruct(chrom, founderAllele);
            //repeat reading pairs of values from hsa and hsb:
            boolean done = false;
            int r = 0;
            while (!done && r<hsbwords.length) {
                if (words[4+r].equals(popdata.missing)) {
                    done = true;
                    if (!hsbwords[r].equals(popdata.missing)) {
                        throw new Exception(fName+".hsa/hsb, line "+(line+1)+": number of non-missing items in both files don't match");
                    }
                }
                else {
                    founderAllele = -1;
                    try {
                        founderAllele = Integer.valueOf(words[r+4]);
                    } catch (Exception ex) {}
                    if (founderAllele<0 || founderAllele>=popdata.founderAlleleCount ) {
                        throw new Exception(fName+".hsa, line "+(line+1)+": invalid founder allele '"+words[r+4]+"'");
                    }
                    double recpos = Double.NaN;
                    try {
                        recpos = Double.valueOf(hsbwords[r]);
                    } catch (Exception ex) {}
                    if (Double.isNaN(recpos)) {
                        throw new Exception(fName+".hsa/hsb, line "+(line+1)+": number of non-missing items in both files don't match");
                    }
                    hast.addSegment(recpos/100.0, founderAllele); //cM to Morgan
                }
                r++;
            }
            //System.out.println("line="+line+"haplostruct:"+hast);
            hs[chrnum][homol] = hast;
            line++;
            if (line % (popdata.chromCount()*popdata.ploidy) == 0) {
                popdata.getIndiv(iName).setHaploStruct(hs);
            }
        } while(true);
        if (line != popdata.indivCount()*popdata.chromCount()*popdata.ploidy) {
            throw new Exception(fName+".hsa/hsb: number of lines incorrect");
        }
    } //readHaploStructFiles

    public static File getParameterFile(String[] args) {
        String filename;
        if (args.length==0) {
            filename = appName+".par";
        } else {
            filename = args[0];
        }
        try {
            File parFile = new File(filename);
            if (parFile.exists() && parFile.isFile()) {
                return parFile;
            }
            else {
                System.out.println(appName+": parameter file "+filename+" not found.");
                return null;
            }
        }
        catch (Exception ex) {
            System.out.println(appName+": parameter file "+filename+" not found.");
            return null;
        }
    } //getParameterFile


    public static String[] getKeyandValue(String line) {
        String[] words = Tools.readWords(line);
        if (words.length==0 || words[0].charAt(0)==commentChar) {
            return new String[] {"",""};
        }
        String key="", value="";
        if (words[0].contains("=")) {
            String[] subwords = words[0].split("=",2);
            key = subwords[0];
            if (subwords.length>1) {
                value = subwords[1];
            }
            else if (words.length>1) {
                value = words[1];
            }
        }
        else {
            key = words[0];
            if (words.length>1) {
                if (words[1].equals("=") && words.length>2) {
                    value = words[2];
                } else {
                    if (words[1].charAt(0)=='=') {
                        value = words[1].substring(1);
                    }
                }
            }
        }
        key = key.toUpperCase();
        return new String[] {key, value};
    } //getKeyandValue



    public static boolean isPositiveInteger(String s) {
        try {
            return Integer.parseInt(s)>0;
        } catch (Exception ex) { return false; }
    }

    public static boolean isInteger(String s) {
        try {
            return Integer.parseInt(s)<=Integer.MAX_VALUE;
        } catch (Exception ex) { return false; }
    }

    public static boolean isProbability(String s) {
        try {
            double d = Double.parseDouble(s);
            return d>=0.0 && d<=1.0;
        } catch (Exception ex) { return false; }
    }

    public static boolean isValidPloidy(String s) {
        try {
            int p = Integer.parseInt(s);
            return p==2 || p==4;
        } catch (Exception ex) { return false; }
    }

    public static boolean isBoolean(String s) {
        s = s.toUpperCase();
        return s.equals("0") || s.equals("1") ||
                s.equals("FALSE") || s.equals("TRUE") ||
                s.equals("NO") || s.equals("YES");
    }

    public static boolean toBoolean(String s) {
        s = s.toUpperCase();
        return s.equals("1") || s.equals("TRUE") || s.equals("YES");
    }

    /**
     * canCreateFile: this will delete an existing file with this name,
     * try to create a new one, and delete it again
     * @param fName
     * @return
     */
    public static boolean canCreateFile(String fName) {
        boolean result = false;
        try {
            File test = new File(fName);
            if (test.exists()) {
                test.delete();
            }
            if (test.createNewFile()) {
                result = test.delete();
            }
        } catch (Exception ex) { }
        return result;
    }

    public static ArrayList<String[]> generateStandardPopulation(String popType,
            int popSize, String missing) {
        ArrayList<String[]> pedList = new ArrayList<String[]>();
        int numlength = Integer.toString(popSize).length();
        DecimalFormat format = new DecimalFormat("0000000000".substring(0, numlength));
        if (popType.equals(strF1)) {
            pedList.add(new String[] {"P1",missing,missing});
            pedList.add(new String[] {"P2",missing,missing});
            for (int i=1; i<=popSize; i++) {
                pedList.add(new String[] {"F1_"+format.format(i),"P1","P2"});
            }
        } else
        if (popType.equals(strF2)) {
            pedList.add(new String[] {"P1",missing,missing});
            pedList.add(new String[] {"P2",missing,missing});
            pedList.add(new String[] {"F1","P1","P2"});
            for (int i=1; i<=popSize; i++) {
                pedList.add(new String[] {"F2_"+format.format(i),"F1","F1"});
            }
        } else
        if (popType.equals(strBC)) {
            pedList.add(new String[] {"P1",missing,missing});
            pedList.add(new String[] {"P2",missing,missing});
            pedList.add(new String[] {"F1","P1","P2"});
            for (int i=1; i<=popSize; i++) {
                pedList.add(new String[] {"BC_"+format.format(i),"P1","F1"});
            }
        } else
        if (popType.equals(strS1)) {
            pedList.add(new String[] {"P1",missing,missing});
            for (int i=1; i<=popSize; i++) {
                pedList.add(new String[] {"S1_"+format.format(i),"P1","P1"});
            }
        }
        // else no population created
        return pedList;
    } //generateStandardPopulation

    public static void simulateHaploStructs(PopulationData popdata) throws Exception {
        for (Individual ind: popdata.getIndividual()) {
            //check is all is well:
            if ( ind.isFounder() && ind.getAllHaploStruct()==null) {
                throw new Exception("Founder was created without haplostruct");
            }
            if ( !ind.isFounder() && ind.getAllHaploStruct()!=null) {
                throw new Exception("Non-founder was created with haplostruct");
            }
            if (!ind.isFounder()) {
                ind.setHaploStruct(ind.getParents()[0].mating(ind.getParents()[1]));
            }
        }
    } //simulateHaploStructs


    public static void testMeiosis(PopulationData popdata, String fName, double gammaFactor) {
        boolean printPubTables = false; //works only with a specific set of chromosomes and maps,
            //therefore not a user-selectable parameter. 
            //Set to false in public distribution
        try {
            
            //preparation:
            popdata.tools.GAMMAFACTOR = gammaFactor;
            PrintWriter out = new PrintWriter(
                    new BufferedWriter(new FileWriter(fName+".dat")));
            out.println("PedigreeSim test output");
            out.println();
            out.println("test settings:");
            String[] settings = popdata.popdataSettings();
            for (int i = 0; i<settings.length; i++) {
                out.println(settings[i]);
            }
            out.println();
            int gamploidy = popdata.ploidy/2;
            String[] locGenotypeNames;
            int[][] locGenotype = null;
            if (gamploidy==1) {
                locGenotypeNames = new String[] {"(0)","(1)"};
            }
            else { //gamploidy==2
                locGenotypeNames = new String[] { "(0,0)","(0,1)","(0,2)","(0,3)",
                    "(1,1)","(1,2)","(1,3)",  "(2,2)","(2,3)", "(3,3)"};
                locGenotype = new int[][] {
                    {0,1,2,3},
                    {1,4,5,6},
                    {2,5,7,8},
                    {3,6,8,9} 
                };
            }
            out.println("chromosome\tlength\tstart\tend\tcentromere\tprefPairing\tfracQuadrivalents");
            for (Chromosome chr : popdata.getChromosome()) {
                out.print(chr.getChromName()+"\t"+(chr.getLength()*100)+"\t"+
                    (chr.getStartPos()*100)+"\t"+(chr.getEndPos()*100)+"\t"+(chr.getCentromerePos()*100));
                int mrkCount = chr.getLocus().size();
                chr.founderAlleleCountAll = new int[popdata.ploidy][mrkCount];
                chr.founderAlleleCountSel = new int[popdata.ploidy][mrkCount];
                chr.locGenotypeCount = new int[locGenotypeNames.length][mrkCount];
                chr.recombCountAll = new int[mrkCount][mrkCount];
                chr.recombCountSel = new int[mrkCount][mrkCount];
                //chr.recombCountCum = new int[mrkCount][mrkCount];
                chr.mapdistCount = new int[mrkCount][mrkCount];
                //chr.recombCountCum = new int[mrkCount][mrkCount];
                chr.recombFrCum = new double[mrkCount][mrkCount];
                chr.recombFrSqCum = new double[mrkCount][mrkCount];
                chr.mapdistCum = new double[mrkCount][mrkCount];
                chr.mapdistSqCum = new double[mrkCount][mrkCount];
                chr.recombPoints = new int[(int)(chr.getLength()*4)];
                chr.founderCount = new int[popdata.ploidy];
                if (popdata.ploidy==4) {
                    TetraploidChromosome chr4 = (TetraploidChromosome)chr;
                    out.print("\t"+chr4.getPrefPairingProb()+"\t"+chr4.getFracQuadrivalents());
                    chr4.founderHomAll = new int[mrkCount];
                    chr4.founderHomSel = new int[mrkCount];
                    chr4.founderHomFrCum = new double[mrkCount];
                    chr4.founderHomFrSqCum = new double[mrkCount];
                    chr4.bivalentCount = 0;
                    chr4.paralQuadrivalentCount = 0;
                    chr4.crossQuadrivalentCount = 0;
                    chr4.exchangeMidFreq = new int[TetraploidChromosome.freqTableLength+1];
                    chr4.exchangeLengthFreq = new int[TetraploidChromosome.freqTableLength+1];
                    chr4.noExchangeLimCount = 0;
                    chr4.oneExchangeLimCount = 0;
                    chr4.twoExchangeLimCount = 0;
                    chr4.exchangeMidSum = 0.0;
                    chr4.exchangeMidSS = 0.0;
                    chr4.exchangeLengthSum = 0.0;
                    chr4.exchangeLengthSS = 0.0;
                    chr4.quadChiasmaSum = 0;
                    chr4.quadChiasmaSS = 0;
                }
                out.println();
            }
            out.flush();
            for (int iter=0; iter<popdata.testIter; iter++) {
                System.out.println("Iteration "+iter);
                //clear the arrays:
                for (Chromosome chr : popdata.getChromosome()) {
                    int mrkCount = chr.getLocus().size();
                    for (int loc=0; loc<mrkCount; loc++) {
                        for (int p=0; p<popdata.ploidy; p++) {
                            chr.founderAlleleCountAll[p][loc]=0;
                            chr.founderAlleleCountSel[p][loc]=0;
                        }
                        for (int loc2=0; loc2<mrkCount; loc2++) {
                            chr.recombCountAll[loc][loc2]=0;
                            chr.recombCountSel[loc][loc2]=0;
                        }
                    }
                    if (popdata.ploidy==4) {
                        for (int loc=0; loc<mrkCount; loc++) {
                            ((TetraploidChromosome)chr).founderHomAll[loc]=0;
                            ((TetraploidChromosome)chr).founderHomSel[loc]=0;
                            //((TetraploidChromosome)chr).founderHomSelSq[loc]=0;
                        }
                    }
                }

                //simulation of gametes:
                if (popdata.testPrintGametes != 0) {
                    out.println("mei\tgam\tchr\thom\tHaploStruct");
                }
                for (int m=0; m<popdata.testMeioseCount; m++) {
                    ArrayList<Gamete> gam = popdata.getIndiv(0).doMeiosis();
                    int selgam = 0; //gamete 0 is a completely random gamete
                    if (popdata.testPrintGametes != 0) {
                        //print all gametes or only selected gamete:
                        for (int g=0; g<gam.size(); g++) {
                            if (g==selgam || popdata.testPrintGametes!=1) {
                                for (int c=0; c<popdata.chromCount(); c++) {
                                    for (int hom=0; hom<popdata.ploidy/2; hom++) {
                                        out.println(m+"\t"+g+"\t"+c+"\t"+hom+"\t"+
                                        gam.get(g).getHaploStruct(c, hom).toString()); 
                                    }    
                                }
                            }    
                        } 
                    }    
                    for (int g=0; g<4; g++) { //loop over 4 gametes in gam
                        for (int c=0; c<popdata.chromCount(); c++) {
                            Chromosome chr = popdata.getChrom(c);
                            int locCount = chr.getLocus().size();
                            //count occurrences of locus genotypes in selected gamete:
                            if (g==selgam) {
                                HaploStruct hs0 = gam.get(g).getHaploStruct(c, 0);
                                for (int loc=0; loc<locCount; loc++) {
                                    double locuspos = chr.getLocus().get(loc).getPosition();
                                    int founderAll0 = hs0.getFounderAt(locuspos);
                                    if (gamploidy==1) {
                                        chr.locGenotypeCount[founderAll0][loc]++; 
                                    }
                                    else { //gamploidy==2
                                        HaploStruct hs1 = gam.get(g).getHaploStruct(c, 1);
                                        int founderAll1 = hs1.getFounderAt(locuspos);
                                        chr.locGenotypeCount[locGenotype[founderAll0][founderAll1]][loc]++;
                                    }            
                                }
                            }

                            //count occurrence of founder alleles and recombinants:
                            for (int p=0; p<gamploidy; p++) {
                                HaploStruct hs = gam.get(g).getHaploStruct(c, p);
                                if (selgam==g) {
                                    chr.founderCount[hs.founderCount()-1] ++;
                                    int i = hs.segmentCount()-1;
                                    if (i>=chr.recombPoints.length) {
                                        int[] temp = new int[chr.recombPoints.length];
                                        System.arraycopy(chr.recombPoints, 0, temp, 0, temp.length);
                                        chr.recombPoints = new int[i+4]; //some extra space
                                        System.arraycopy(temp, 0, chr.recombPoints, 0, temp.length);
                                    }
                                    chr.recombPoints[i] ++;
                                }    
                                for (int loc=0; loc<locCount; loc++) {
                                    double locuspos = chr.getLocus().get(loc).getPosition();
                                    int locFounderAllele = hs.getFounderAt(locuspos);
                                    chr.founderAlleleCountAll[locFounderAllele][loc]++;
                                    if (selgam==g) {
                                        chr.founderAlleleCountSel[locFounderAllele][loc]++;
                                    }    
                                    for (int loc2=0; loc2<loc; loc2++) {
                                        int loc2FounderAllele =
                                            hs.getFounderAt(chr.getLocus().get(loc2).getPosition());
                                        if (loc2FounderAllele != locFounderAllele) {
                                            chr.recombCountAll[loc][loc2]++;
                                            if (selgam==g) {
                                                chr.recombCountSel[loc][loc2]++;
                                            }
                                        }
                                    }
                                }
                            } //for p

                            //count occurrences of homozygosity:
                            if (popdata.ploidy==4) {
                                boolean[] homoz = gam.get(g).homozygous(c);
                                for (int loc=0; loc<locCount; loc++) {
                                    ((TetraploidChromosome)(chr)).founderHomAll[loc]
                                            += homoz[loc] ? 1 : 0;
                                    if (selgam==g) {
                                        ((TetraploidChromosome)(chr)).founderHomSel[loc]
                                            += homoz[loc] ? 1 : 0;
                                    }
                                }
                            }
                        } // for c
                    } //for g
                } // for m

                //output of cumulative results for current iteration:
                if (popdata.testPrintEachIter) {
                    if (popdata.testIter>1) {
                        out.println();
                        out.println("Iteration "+iter);
                    }
                    out.println();
                    out.println("Cumulative results for ALL gametes per meiosis:");
                    double total = 4 * gamploidy * popdata.testMeioseCount;

                    for (Chromosome chr: popdata.getChromosome()) {
                        out.println();
                        out.println("Chromosome "+chr.getChromName());
                        
                        out.println();
                        out.println("Founder allele frequencies:");
                        chr.printMapHorizontal(out,true);
                        for (int p=0; p<popdata.ploidy; p++) {
                            out.print(p);
                            for (int loc=0; loc<chr.getLocus().size(); loc++)
                                out.print("\t"+(chr.founderAlleleCountAll[p][loc]/total));
                            out.println();
                        }

                        out.println();
                        out.println("Recombination frequencies:");
                        for (int loc=0; loc<chr.getLocus().size(); loc++) {
                            out.print(chr.getLocus().get(loc).getLocusName());
                            for (int loc2=0; loc2<loc; loc2++)
                                out.print("\t"+(chr.recombCountAll[loc][loc2]/total));
                            out.println();
                        }
                        chr.printMapHorizontal(out, false);
                        
                        if (popdata.testPrintMapdistances) {    
                            out.println();
                            out.println("Inferred map distances:");
                            for (int loc=0; loc<chr.getLocus().size(); loc++) {
                                out.print(chr.getLocus().get(loc).getLocusName());
                                for (int loc2=0; loc2<loc; loc2++)
                                    out.print("\t"+(Tools.recombToMapdist(chr.recombCountAll[loc][loc2]/total,
                                            popdata.chiasmaInterference)));
                                out.println();
                            }
                            chr.printMapHorizontal(out, false);
                        }
                        
                        if (popdata.ploidy==4) {
                            total = 4 * popdata.testMeioseCount;
                            out.println();
                            out.println("Homozygosity frequencies for founder alleles:");
                            chr.printMapHorizontal(out, true);
                            out.print("freq");
                            for (int loc=0; loc<chr.getLocus().size(); loc++)
                                out.print("\t"+((double)(((TetraploidChromosome)(chr)).founderHomAll[loc])/total));
                            out.println();
                        }
                    } //for chr

                    out.println();
                    out.println("Cumulative results for ONE RANDOM gamete per meiosis:");
                    total = gamploidy * popdata.testMeioseCount;

                    for (Chromosome chr: popdata.getChromosome()) {
                        out.println();
                        out.println("Chromosome "+chr.getChromName());
                        out.println();
                        out.println("Founder allele frequencies:");
                        chr.printMapHorizontal(out, true);
                        for (int p=0; p<popdata.ploidy; p++) {
                            out.print(p);
                            for (int loc=0; loc<chr.getLocus().size(); loc++)
                                out.print("\t"+(chr.founderAlleleCountSel[p][loc]/total));
                            out.println();
                        }
                        
                        out.println();
                        out.println("Recombination frequencies:");
                        for (int loc=0; loc<chr.getLocus().size(); loc++) {
                            out.print(chr.getLocus().get(loc).getLocusName());
                            for (int loc2=0; loc2<loc; loc2++)
                                out.print("\t"+(chr.recombCountSel[loc][loc2]/total));
                            out.println();
                        }
                        chr.printMapHorizontal(out, false);
                        
                        if (popdata.testPrintMapdistances) {
                            out.println();
                            out.println("Inferred map distances:");
                            double recomb, dist;
                            for (int loc=0; loc<chr.getLocus().size(); loc++) {
                                out.print(chr.getLocus().get(loc).getLocusName());
                                for (int loc2=0; loc2<loc; loc2++) {
                                    recomb = chr.recombCountSel[loc][loc2]/total;
                                    dist = Tools.recombToMapdist(recomb,popdata.chiasmaInterference);
                                    out.print("\t"+dist);
                                }
                                out.println();
                            }
                            chr.printMapHorizontal(out, false);
                        }

                        if (popdata.ploidy==4) {
                            total = popdata.testMeioseCount;
                            out.println();
                            out.println("Homozygosity frequencies for founder alleles:");
                            chr.printMapHorizontal(out, true);
                            out.print("freq");
                            for (int loc=0; loc<chr.getLocus().size(); loc++)
                                out.print("\t"+((double)(((TetraploidChromosome)(chr)).founderHomSel[loc])/total));
                            out.println();
                        }
                    } //for chr
                } //only output per iter if popdata.testIter<=maxOut

                //add the results for this iter (for the one selected gamete per meiosis) to the totals:
                double total = gamploidy * popdata.testMeioseCount;
                for (Chromosome chr: popdata.getChromosome()) {
                    //add the recomb and map distance results for this iter to the totals
                    double recomb, dist;
                    for (int loc=0; loc<chr.getLocus().size(); loc++) {
                        for (int loc2=0; loc2<loc; loc2++) {
                            recomb = chr.recombCountSel[loc][loc2]/total;
                            chr.recombFrCum[loc][loc2] += recomb;
                            chr.recombFrSqCum[loc][loc2] += recomb*recomb;
                            dist = Tools.recombToMapdist(recomb,popdata.chiasmaInterference); //not correct for quadrivalents!
                            if (!Double.isNaN(dist)) {
                                chr.mapdistCount[loc][loc2] ++;
                                chr.mapdistCum[loc][loc2] += dist;
                                chr.mapdistSqCum[loc][loc2] += dist*dist;
                            }
                        }
                    }
                    //add the homozygosity frequency for this iter to the totals
                    if (popdata.ploidy==4) {
                        for (int i=0; i<((TetraploidChromosome)(chr)).founderHomSel.length; i++) {
                            double homozFreq =
                                ((double)((TetraploidChromosome)(chr)).founderHomSel[i])
                                    / popdata.testMeioseCount;
                            ((TetraploidChromosome)(chr)).founderHomFrCum[i] +=
                                      homozFreq;
                            ((TetraploidChromosome)(chr)).founderHomFrSqCum[i] +=
                                      homozFreq*homozFreq;
                        }
                    }
                } //for chr
            } //for iter

            //output accumulated recombination statistics over all iterations:
            int recombfunction=0; //calculated per chromosome; for publication tables we assume
                      //that all chromosomes have the same recombfunction
            TDistribution tDistribution = new TDistribution(popdata.testIter-1);
            out.println();
            out.println("Accumulated statistics over all "+popdata.testIter+" iterations:");
            for (Chromosome chr: popdata.getChromosome()) {
                out.println();
                out.println("Chromosome "+chr.getChromName());

                out.println();
                out.println("Frequency distribution of the number of recombination points");
                int mi = chr.recombPoints.length-1;
                while (mi>=0 && chr.recombPoints[mi]==0) mi--;
                int sum=0, count=0;
                out.println("points\tfrequency");
                for (int i=0; i<=mi; i++) {
                    out.println(i+"\t"+chr.recombPoints[i]);
                    count += chr.recombPoints[i];
                    sum += i*chr.recombPoints[i];
                }
                out.println("count & mean\t"+count+"\t"+(((double)sum)/count));

                out.println();
                out.println("Frequency distribution of the number of different founder alleles");
                out.println("alleles\tfrequency");
                sum=0; count=0;
                for (int i=0; i<popdata.ploidy; i++) {
                    out.println((i+1)+"\t"+chr.founderCount[i]);
                    count += chr.founderCount[i];
                    sum += (i+1)*chr.founderCount[i];
                }
                out.println("count & mean\t"+count+"\t"+(((double)sum)/count));
                
                out.println();
                out.println("Frequency distribution of locus genotypes");
                chr.printMapHorizontal(out, true);
                for (int i=0; i<locGenotypeNames.length; i++) {
                    out.print(locGenotypeNames[i]);
                    for (int loc=0; loc<chr.getLocus().size(); loc++) {
                        out.print("\t"+chr.locGenotypeCount[i][loc]);
                    }
                    out.println();
                }

                int mainconfig = 0; //0: bivalents
                int mainconfigCount = popdata.testIter*popdata.testMeioseCount; //for diploids
                if (popdata.ploidy==4) {
                    out.println();
                    out.println("Homozygosity frequencies for founder alleles:");
                    chr.printMapHorizontal(out, true);
                    out.print("freq");
                    for (int loc=0; loc<chr.getLocus().size(); loc++)
                        out.print("\t"+(((TetraploidChromosome)(chr)).founderHomFrCum[loc])
                                / popdata.testIter);
                    out.println();
                    out.print("stdev");
                    for (int loc=0; loc<chr.getLocus().size(); loc++)
                        out.print("\t"+Tools.sampleStDev(
                                ((TetraploidChromosome)(chr)).founderHomFrCum[loc], 
                                ((TetraploidChromosome)(chr)).founderHomFrSqCum[loc], 
                                popdata.testIter));
                    out.println();

                    out.println();
                    out.println("Total number of two-bivalent configs\t" +
                            ((TetraploidChromosome)(chr)).bivalentCount/2);
                    out.println("Total number of parallel quadrivalents\t" +
                            ((TetraploidChromosome)(chr)).paralQuadrivalentCount);
                    out.println("Total number of cross-type quadrivalents\t" +
                            ((TetraploidChromosome)(chr)).crossQuadrivalentCount);
                    mainconfig = 0; //0: 2 bivalents
                    mainconfigCount = ((TetraploidChromosome)(chr)).bivalentCount/2;
                    if (((TetraploidChromosome)(chr)).paralQuadrivalentCount >
                            mainconfigCount) {
                        mainconfig = 1; //1: parallel quadrivalents
                        mainconfigCount = ((TetraploidChromosome)(chr)).paralQuadrivalentCount;
                    }
                    if (((TetraploidChromosome)(chr)).crossQuadrivalentCount >
                            mainconfigCount) {
                        mainconfig = 2; //2: cross-type quadrivalents
                        mainconfigCount = ((TetraploidChromosome)(chr)).crossQuadrivalentCount;
                    } 
                    out.println("Data on the chromosome exchange interval of cross-type quadrivalents:");
                    out.println("Number of cross-type quadrivalents where the chromosome exchange interval is:");
                    out.println("not defined (no chiasmata):\t"+
                            ((TetraploidChromosome)(chr)).noExchangeLimCount);
                    out.println("defined on one side only:\t" +
                            ((TetraploidChromosome)(chr)).oneExchangeLimCount);
                    out.println("defined on both sides:\t" +
                            ((TetraploidChromosome)(chr)).twoExchangeLimCount);
                    out.println("Distribution of midpoints of cross-type quadrivalent exchange intervals (only if chiasmata on both sides)");
                    out.println("interval (M)\t\t\tfrequency");
                    for (int i=0; i<((TetraploidChromosome)(chr)).exchangeMidFreq.length; i++) {
                        out.println((chr.getStartPos()+i*(chr.getLength()/TetraploidChromosome.freqTableLength)) +
                                "\t<=X<\t" +
                                (chr.getStartPos()+(i+1)*(chr.getLength()/TetraploidChromosome.freqTableLength)) +
                                "\t" + (((double)(((TetraploidChromosome)(chr)).exchangeMidFreq[i]))/
                                 ((TetraploidChromosome)(chr)).twoExchangeLimCount));
                    }
                    double mean =  ((TetraploidChromosome)(chr)).exchangeMidSum /
                            ((TetraploidChromosome)(chr)).twoExchangeLimCount;
                    out.println("mean\t" + mean);
                    out.println("st.dev.\t" + Tools.sampleStDev(
                            ((TetraploidChromosome)(chr)).exchangeMidSum,
                            ((TetraploidChromosome)(chr)).exchangeMidSS,
                            ((TetraploidChromosome)(chr)).twoExchangeLimCount));
                    out.println("Distribution of lengths of cross-type quadrivalent exchange intervals (only if chiasmata on both sides)");
                    out.println("interval\t\t\tfrequency");
                    //Note: with popdata.quadriMethod==2 and no chiasma interference the length is always 0: 
                    //the last failed chiasma extends the exchangeLim to the opposing chiasma or chromosome end
                    for (int i=0; i<((TetraploidChromosome)(chr)).exchangeLengthFreq.length; i++) {
                        out.println((chr.getStartPos()+i*(chr.getLength()/TetraploidChromosome.freqTableLength)) +
                                "\t<=X<\t" +
                                (chr.getStartPos()+(i+1)*(chr.getLength()/TetraploidChromosome.freqTableLength)) +
                                "\t" + (((double)(((TetraploidChromosome)(chr)).exchangeLengthFreq[i]))/
                                 ((TetraploidChromosome)(chr)).twoExchangeLimCount));
                    }
                    mean =  ((TetraploidChromosome)(chr)).exchangeLengthSum /
                            ((TetraploidChromosome)(chr)).twoExchangeLimCount;
                    out.println("mean\t" + mean);
                    out.println("st.dev.\t" + Tools.sampleStDev(
                            ((TetraploidChromosome)(chr)).exchangeLengthSum,
                            ((TetraploidChromosome)(chr)).exchangeLengthSS,
                            ((TetraploidChromosome)(chr)).twoExchangeLimCount));
                } //if ploidy==4

                out.println();
                out.print("Expected recombination frequencies: ");
                recombfunction=0; //bivalents
                if (popdata.ploidy==2 || 
                    mainconfigCount<popdata.testIter*popdata.testMeioseCount ||
                    mainconfig==0) {
                    if (popdata.chiasmaInterference) out.println("Kosambi");
                    else out.println("Haldane");
                }
                else if (mainconfig==1) {
                    recombfunction=1;
                    if (popdata.chiasmaInterference) out.println("Sved_parallel based on Kosambi");
                    else out.println("Sved_parallel based on Haldane");
                }
                else {
                    assert mainconfig==2;
                    recombfunction=2;
                    if (popdata.chiasmaInterference) out.println("Sved_cross-type based on Kosambi");
                    else out.println("Sved_cross-type based on Haldane");
                }
                for (int loc=0; loc<chr.getLocus().size(); loc++) {
                    out.print(chr.getLocus().get(loc).getLocusName());
                    for (int loc2=0; loc2<loc; loc2++) {
                        double dist = chr.getLocus().get(loc).getPosition() -
                                chr.getLocus().get(loc2).getPosition();
                        double expRec = recombfunction==0 ? 
                                Tools.mapdistToRecomb(dist,popdata.chiasmaInterference) :
                                recombfunction==1 ?
                                Tools.recombBivalentToRecombParallelQuadrivalent(
                                    Tools.mapdistToRecomb(dist,popdata.chiasmaInterference)) :
                                Tools.recombBivalentToRecombCrossQuadrivalent(
                                    Tools.mapdistToRecomb(dist,popdata.chiasmaInterference));
                        out.print("\t"+expRec);
                    }
                    out.println();
                }
                chr.printMapHorizontal(out, false);

                out.println();
                out.println("Observed Recombination frequencies:");
                for (int loc=0; loc<chr.getLocus().size(); loc++) {
                    out.print(chr.getLocus().get(loc).getLocusName());
                    for (int loc2=0; loc2<loc; loc2++) {
                        double rec = chr.recombFrCum[loc][loc2]/popdata.testIter;
                        out.print("\t"+rec);
                    }
                    out.println();
                }
                chr.printMapHorizontal(out, false);

                out.println();
                out.println("Deviations from expected recombination frequencies:");
                for (int loc=0; loc<chr.getLocus().size(); loc++) {
                    out.print(chr.getLocus().get(loc).getLocusName());
                    for (int loc2=0; loc2<loc; loc2++) {
                        double rec = chr.recombFrCum[loc][loc2]/popdata.testIter;
                        double dist = chr.getLocus().get(loc).getPosition() -
                                chr.getLocus().get(loc2).getPosition();
                        double expRec = recombfunction==0 ? 
                                Tools.mapdistToRecomb(dist,popdata.chiasmaInterference) :
                                recombfunction==1 ?
                                Tools.recombBivalentToRecombParallelQuadrivalent(
                                    Tools.mapdistToRecomb(dist,popdata.chiasmaInterference)) :
                                Tools.recombBivalentToRecombCrossQuadrivalent(
                                    Tools.mapdistToRecomb(dist,popdata.chiasmaInterference));
                        out.print("\t"+(rec-expRec));
                    }
                    out.println();
                }
                chr.printMapHorizontal(out, false);

                out.println();
                out.println("Standard Deviations (sample) of recombination frequencies:");
                for (int loc=0; loc<chr.getLocus().size(); loc++) {
                    out.print(chr.getLocus().get(loc).getLocusName());
                    for (int loc2=0; loc2<loc; loc2++) {
                        out.print("\t"+Tools.sampleStDev(
                                chr.recombFrCum[loc][loc2],
                                chr.recombFrSqCum[loc][loc2],
                                popdata.testIter));
                    }
                    out.println();
                }
                chr.printMapHorizontal(out, false);

                out.println();
                out.println("Significance of deviations from expected recombination frequencies:");
                double stdev, tvalue, tprob;
                for (int loc=0; loc<chr.getLocus().size(); loc++) {
                    out.print(chr.getLocus().get(loc).getLocusName());
                    for (int loc2=0; loc2<loc; loc2++) {
                        double rec = chr.recombFrCum[loc][loc2]/popdata.testIter;
                        double dist = chr.getLocus().get(loc).getPosition() -
                                chr.getLocus().get(loc2).getPosition();
                        double expRec = recombfunction==0 ? 
                                Tools.mapdistToRecomb(dist,popdata.chiasmaInterference) :
                                recombfunction==1 ?
                                Tools.recombBivalentToRecombParallelQuadrivalent(
                                    Tools.mapdistToRecomb(dist,popdata.chiasmaInterference)) :
                                Tools.recombBivalentToRecombCrossQuadrivalent(
                                    Tools.mapdistToRecomb(dist,popdata.chiasmaInterference));
                        stdev= Tools.sampleStDev(
                                chr.recombFrCum[loc][loc2],
                                chr.recombFrSqCum[loc][loc2],
                                popdata.testIter);
                        tvalue = -Math.abs(rec-expRec)/(stdev/Math.sqrt(popdata.testIter));
                        tprob = tDistribution.cumulative(tvalue)*2.0; //two-sided test against expectation
                        out.print("\t"+tprob);
                    }
                    out.println();
                }
                chr.printMapHorizontal(out, false);

                if (popdata.testPrintMapdistances) {
                    out.println();
                    out.println("number of iterations with recombination >= 0.5 :");
                    for (int loc=0; loc<chr.getLocus().size(); loc++) {
                        out.print(chr.getLocus().get(loc).getLocusName());
                        for (int loc2=0; loc2<loc; loc2++) {
                            out.print("\t"+(popdata.testIter-chr.mapdistCount[loc][loc2]));
                        }
                        out.println();
                    }
                    chr.printMapHorizontal(out, false);

                    out.println();
                    out.println("Average inferred map distances (based on bivalent map function):");
                    out.println();
                    for (int loc=0; loc<chr.getLocus().size(); loc++) {
                        out.print(chr.getLocus().get(loc).getLocusName());
                        for (int loc2=0; loc2<loc; loc2++) {
                            double meandist = chr.mapdistCum[loc][loc2]/chr.mapdistCount[loc][loc2];
                            out.print("\t"+meandist);
                        }
                        out.println();
                    }
                    chr.printMapHorizontal(out, false);

                    out.println();
                    out.println("Deviations from true distances (based on bivalent map function):");
                    for (int loc=0; loc<chr.getLocus().size(); loc++) {
                        out.print(chr.getLocus().get(loc).getLocusName());
                        for (int loc2=0; loc2<loc; loc2++) {
                            double meandist = chr.mapdistCum[loc][loc2]/chr.mapdistCount[loc][loc2];
                            double expDist = chr.getLocus().get(loc).getPosition() -
                                    chr.getLocus().get(loc2).getPosition();
                            out.print("\t"+(meandist-expDist));
                        }
                        out.println();
                    }
                    chr.printMapHorizontal(out, false);

                    out.println();
                    out.println("Standard Deviations (sample) of map distances:");
                    for (int loc=0; loc<chr.getLocus().size(); loc++) {
                        out.print(chr.getLocus().get(loc).getLocusName());
                        for (int loc2=0; loc2<loc; loc2++) {
                            out.print("\t"+Tools.sampleStDev(
                                    chr.mapdistCum[loc][loc2],
                                    chr.mapdistSqCum[loc][loc2],
                                    chr.mapdistCount[loc][loc2]));
                        }
                        out.println();
                    }
                    chr.printMapHorizontal(out, false);

                    out.println();
                    out.println("Significance of deviations from true distances (based on bivalent map function):");
                    for (int loc=0; loc<chr.getLocus().size(); loc++) {
                        out.print(chr.getLocus().get(loc).getLocusName());
                        for (int loc2=0; loc2<loc; loc2++) {
                            double meandist = chr.mapdistCum[loc][loc2]/chr.mapdistCount[loc][loc2];
                            stdev= Tools.sampleStDev(
                                    chr.mapdistCum[loc][loc2],
                                    chr.mapdistSqCum[loc][loc2],
                                    chr.mapdistCount[loc][loc2]);
                            double expDist = chr.getLocus().get(loc).getPosition() -
                                    chr.getLocus().get(loc2).getPosition();
                            tvalue = -Math.abs(meandist-expDist)/(stdev/Math.sqrt(chr.mapdistCount[loc][loc2]));
                            if (chr.mapdistCount[loc][loc2]==popdata.testIter)
                                tprob = tDistribution.cumulative(tvalue)*2.0; //two-sided test against expectation
                            else if (chr.mapdistCount[loc][loc2]<2)
                                tprob = Double.NaN;
                            else tprob = 1.0 - new TDistribution(chr.mapdistCount[loc][loc2]).cumulative(tvalue);
                            out.print("\t"+tprob);
                        }
                        out.println();
                    }
                    chr.printMapHorizontal(out, false);
                }

            } // for chr
            out.flush();
            
            if (printPubTables) {
                //Note: this only works with a specific set of
                //chromosomes and maps
                out.println();
                out.print("Recombination table; expected recombination based on ");
                if (recombfunction==0) {
                    if (popdata.chiasmaInterference) out.println("Kosambi");
                    else out.println("Haldane");
                }
                else if (recombfunction==1) {
                    if (popdata.chiasmaInterference) out.println("Sved_parallel based on Kosambi");
                    else out.println("Sved_parallel based on Haldane");
                }
                else {
                    assert recombfunction==2;
                    if (popdata.chiasmaInterference) out.println("Sved_cross-type based on Kosambi");
                    else out.println("Sved_cross-type based on Haldane");
                }    
                int[][][][] recloc = new int[][][][] 
                {
                    { //chrom A
                        { {0,1},{0,2},{0,3},{0,4},{0,6} },
                        { {9,11},{8,12},{7,13},{6,14},{5,15},{0,20} },
                        { {19,20},{18,20},{17,20},{16,20},{14,20} } 
                    }, 
                    { //chrom B
                        { {0,1},{0,2},{0,3},{0,4},{0,5},{0,12} },
                        { {11,13},{10,14},{9,15},{8,16},{7,17},{6,18},{0,24} },
                        { {23,24},{22,24},{21,24},{20,24},{19,24},{12,24} }
                    },
                    { //chrom C
                        { {0,1},{0,2},{0,3},{0,4},{0,5},{0,7},{0,15} },
                        { {14,16},{13,17},{12,18},{11,19},{10,20},{8,22},{7,23},{0,30} },
                        { {29,30},{28,30},{27,30},{26,30},{25,30},{23,30},{15,30} }
                    },
                    { //chrom D
                        { {0,1},{0,2},{0,3},{0,4},{0,5},{0,6},{0,8},{0,17},{0,34} },
                        { {16,18},{15,19},{14,20},{13,21},{12,22},{11,23},{9,25},{8,26},{0,34} },
                        { {33,34},{32,34},{31,34},{30,34},{29,34},{28,34},{26,34},{17,34} }
                    }
                };
                out.println("\t\tA\t\tA\t\tA\t\tB\t\tB\t\tB\t\tC\t\tC\t\tC\t\tD\t\tD\t\tD");
                out.println("interval\trec\ttop\t\tcenter\t\tbottom\t\ttop\t\tcenter\t\tbottom\t\ttop\t\tcenter\t\tbottom\t\ttop\t\tcenter\t\tbottom");
                for (int line=0; line<9; line++) {
                    double dist = popdata.getChrom(3).getLocus().get(recloc[3][0][line][1]).position;
                    double expRec = recombfunction==0 ? 
                            Tools.mapdistToRecomb(dist,popdata.chiasmaInterference) :
                            recombfunction==1 ?
                            Tools.recombBivalentToRecombParallelQuadrivalent(
                                Tools.mapdistToRecomb(dist,popdata.chiasmaInterference)) :
                            Tools.recombBivalentToRecombCrossQuadrivalent(
                                Tools.mapdistToRecomb(dist,popdata.chiasmaInterference));
                    //double expRec = Tools.mapdistToRecomb(dist,popdata.chiasmaInterference);
                    out.print(dist+"\t"+expRec);
                    for (int c=0; c<4; c++) {
                        Chromosome ch = popdata.getChrom(c);
                        for (int col=0; col<3; col++) {
                            if (recloc[c][col].length>line) {
                                int locA = recloc[c][col][line][0];
                                int locB = recloc[c][col][line][1];
                                double rec = ch.recombFrCum[locB][locA]/popdata.testIter;
                                double stdev= Tools.sampleStDev(
                                        ch.recombFrCum[locB][locA],
                                        ch.recombFrSqCum[locB][locA],
                                        popdata.testIter);
                                double tvalue = -Math.abs(rec-expRec)/(stdev/Math.sqrt(popdata.testIter));
                                double tprob = tDistribution.cumulative(tvalue)*2.0; //two-sided test against expectation
                                out.print("\t"+rec+"\t"+tprob);
                            }
                            else out.print("\t\t");
                        }    
                    }
                    out.println();
                } //for line
            } //printPubTables

            out.flush();
        } catch (Exception ex) {
            System.out.println(ex.getMessage());
        }
    } //testMeiosis

}