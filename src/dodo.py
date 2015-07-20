#########################################################################
# Jim Kaminski
# Huttenhower Lab
# ShortBRED Paper DoIT Script
# 7/3/2015
#########################################################################

"""
This is a Python doit script which creates the synthetoc metagenomes for the
2015 ShortBRED paper. Originally, this was performed by sfle in the
repositories shortbred-mg and shortbred-sfle, but was recoded into this
project to make running the analysis easier.

Dependencies Required:
Python 2.7
Python packages: doit, Biopython, numpy, scipy
R Packages: optparse,

Additional Actions:
* Please set "dirName" to the directory where you downloaded this project.
* Download ShortBRED, place inside "src".
* Please note you will need to edit the alias file in the blastdb,
"refdb.pal", to point to its location on your system.
* If running in hutlab, use the following command at the shell:
export R_LIBS=$R_LIBS:"/n/home11/jkaminski/Rlibs"


"""
import os, glob, sys, re
import math

#############################################################################
# Functions

def python_check_create_dir( strDir ):
    if not os.path.exists( strDir ):
        os.makedirs( strDir )
    return strDir

#############################################################################
# Constants


dirMain = "/n/huttenhower_lab/data/shortbred/shortbred_doit"
c_blastdbRefName = "/n/huttenhower_lab/data/shortbred/blastdb/refdb"


c_cmdBlastP			=   "/n/home11/jkaminski/blast/ncbi-blast-2.2.28+/bin/blastp"
c_cmdMuscle          = "/n/huttenhower_lab/tools/muscle/muscle3.8.31_i86linux64"
c_cmdCDhit           = "/n/huttenhower_lab/tools/cdhit/cd-hit"
c_cmdUsearch= "/n/home11/jkaminski/usearch/usearch"
c_cmdMakeBlast = "/n/home11/jkaminski/blast/ncbi-blast-2.2.28+/bin/makeblastdb"

c_dirBigTmp =  "/n/hutlab12_nobackup/data"


# sandy version (Added to recreate settings for when we ran this on Sandy in
# 2013 and 2014.)
c_cmdSandyBlastP			=    "/n/home11/jkaminski/blast/ncbi-blast-2.2.28+/bin/blastp"
c_cmdSandyMuscle          = "/n/sw/muscle-3.8.31/bin/muscle"
c_cmdSandyCDhit           = "/n/sw/cd-hit-v4.6-2012-04-25/cd-hit"
c_cmdUsearch= "/n/home11/jkaminski/usearch/usearch"
c_cmdSandyMakeBlast = "/n/sw/ncbi-blast-2.2.28/bin/makeblastdb"



dirData =  dirMain + os.sep + "data"
dirOut =  dirMain + os.sep + "out"
dirTmp =  dirMain + os.sep + "tmp"
dirSrc = dirMain + os.sep + "src"
dirFigures = dirOut + os.sep + "figures"

dirShortBRED = dirSrc + os.sep + "shortbred"
dirBacteriaGenomes =  dirData + os.sep + "bacteria_genomes"
dirMGSpecs = dirData + os.sep + "synthetic_specs"
txtGenomes = dirBacteriaGenomes + os.sep + "GenomeNames.txt"

pyCheckSeqs = dirShortBRED + os.sep + "src" + os.sep + "check_sequences.py"
pyMutateNucs = dirSrc + os.sep + "create_metagenome" + os.sep + "MutateNucs.py"
pyRemovePosGenes = dirSrc + os.sep + "create_metagenome" + os.sep + "RemovePosGenes.py"
pySimpleSim = dirSrc + os.sep + "create_metagenome" + os.sep + "simplesim.py"
pyGemSim = dirSrc + os.sep + "create_metagenome" + os.sep + "gemsim" + os.sep + "GemReads.py"
pyCheckLengths = dirSrc + os.sep + "create_metagenome" + os.sep + "check_seq_length.py"
pyFastq2Fasta = dirSrc + os.sep + "create_metagenome" + os.sep + "fastq2fasta.py"
pyCountWGS = dirSrc + os.sep + "utils" + os.sep + "getWGScounts.py"
pyCheckGenes = dirSrc + os.sep + "utils" + os.sep + "checkgene.py"
pyScreenNucs = dirSrc + os.sep + "utils" + os.sep + "ScreenNucs.py"

pySBIdentify = dirShortBRED + os.sep + "shortbred_identify.py"
pySBQuantify = dirShortBRED + os.sep + "shortbred_quantify.py"

rEval = dirSrc + os.sep + "utils" + os.sep + "REval.R"
rMergeData = dirSrc + os.sep + "utils" + os.sep + "MergeFiles.R"
rPlotAllAUC = dirSrc + os.sep + "utils" + os.sep + "AUCPlot.R"
rFig2 = dirSrc + os.sep + "utils" + os.sep + "MakeFigure2.R"

c_iBlastThreads = 20
c_dConsThresh = 0.95
c_dMinIDForMatch = .90
c_dMatchMaxLen = .15
c_dQuasiClust = .85
c_iQuasiThresh = 1
c_dMatchID = .95
c_iCentLen = 30
c_iTotalMarkerLength = 300
c_iXLimit = 1




strSBIdentifyParams =  " --consthresh " +str(c_dConsThresh) + " --blastp " + c_cmdBlastP + " --muscle " + c_cmdMuscle  + " --cdhit " + c_cmdCDhit + " --usearch " + c_cmdUsearch + " --makeblastdb " + c_cmdMakeBlast


astrDB = ["ARDB","VF"]
"""
###################################
Create Synthetic Metagenomes
###################################
* Create Nucleotide Sequences Based on Input Proteins (may mutate)
* Screen Bacterial Genomes and Remove Matches to Input Sequences
* Make Abundance Table and Gold Standard, add material to ends of nucleotide sequences
* Make synthetic metagenomes with GemSIM (fastq)
* Convert to fasta
* Reduce gene lables for spiked genes to ">USR_NAME_END"
* Count reads
"""
c_dMutRate = 0.03
dirSynthMG = dirOut + os.sep + "synthetic_metagenomes"
dirCleanGenomes = dirSynthMG + os.sep + "clean_genomes"
python_check_create_dir(dirSynthMG)
python_check_create_dir(dirCleanGenomes)
astrGenomes =  ["b.thetaiotaomicron.genome", "p.gingivalis.genome", "p.distasonis.genome",
"p.ruminicola.genome", "v.parvula.genome", "b.mallei.genome", "l.acidophilus.genome",
"e.eligens.genome", "s.aureus_n315.genome", "e.coli.genome", "b.longum.genome","f.nucleatum.genome",
"s.pneumoniae.genome", "c.difficile.genome", "n.meningitidis.genome", "c.jejuni.genome","p.multocida.genome",
"f.johnsoniae.genome","l.buccalis.genome","r.mucilaginosa.genome"]


astrMG = ["Illumina","Illumina-mutated","Illumina-mutated-5pct","454_mockYAT","454_deep","small"]
astrPct = ["05","10","25"]

def task_PrepareNucs():
    for strDB in ["ARDB","VF"]:

                python_check_create_dir(dirCleanGenomes + os.sep + strDB)
                faaInput = dirData + os.sep + "input_sequences" + os.sep + strDB + ".input.faa"
                faaCleanProteinSeqs = dirSynthMG + os.sep +strDB + "clean.faa"
                fnaMutNucs = dirSynthMG + os.sep +strDB + "mutated.fna"


                yield {'name': "Screen input prot file for nucs: " +strDB ,
                   'actions': ["python " + pyCheckSeqs + " < " + faaInput + " > " + faaCleanProteinSeqs],
                   'targets': [faaCleanProteinSeqs],
                   'file_dep': [faaInput],}



                for strMG in astrMG:


                    if strMG == "Illumina":
                        iReads = 5000000
                        iReadLen = "d" # "d" tells GemSim ("GemReads.py") to use the empirical model for length
                        zipModel = dirMGSpecs + os.sep + "ill100v5_s.gzip" # Illumina Model
                        iPadLength = 100
                    elif strMG== "454_mockYAT":
                        iReads = 155890
                        iReadLen = "d"
                        zipModel = dirMGSpecs + os.sep + "r454ti_s.gzip" # 454 Model
                        iPadLength = 450
                    elif strMG=="454_deep":
                        iReads = 1250000
                        iReadLen = "d"
                        zipModel = dirMGSpecs + os.sep + "r454ti_s.gzip"  # 454 Model
                        iPadLength = 450
                    elif strMG=="small":
                        iReads = 4000000
                        iReadLen = "d"
                        zipModel = dirMGSpecs + os.sep + "ill100v5_s.gzip"  # 454 Model
                        iPadLength = 100

                    for strRun in astrPct:
                        if (strRun == "05"):
                            iGenes = 150
                            dSpike = int(strRun)/100.0
                        if (strRun == "10"):
                            iGenes = 500
                            dSpike = int(strRun)/100.0
                        elif (strRun == "25"):
                            iGenes = 1000
                            dSpike = int(strRun)/100.0






                        dirOutMG = dirSynthMG + os.sep + strDB + os.sep + strMG + "-" + strRun
                        dirFinalGenes = dirOutMG + os.sep + "final_genes"
                        python_check_create_dir(dirSynthMG  + os.sep + strDB)
                        python_check_create_dir(dirSynthMG  + os.sep + strDB + os.sep + strMG)
                        python_check_create_dir(dirOutMG)
                        python_check_create_dir(dirFinalGenes)
                        fastqSim = dirOutMG +os.sep + strMG + "-" + strRun+ "_single.fastq"
                        fastaSim = dirOutMG +os.sep + strMG + "-" + strRun+ ".fasta"

                        # Choose appropriate nucelotide spike in.
                        if strMG == "Illumina-mutated":
                            fnaMutNucs = dirOutMG + os.sep +strDB + "mutated.fna"
                            yield {'name': "Screen input prot file for nucs, make mutated nucs: " +strDB + strRun + strMG,
                            'actions': ["python " + pyCheckSeqs + " < " + faaInput + " | " + pyMutateNucs + " --mutrate " + str(c_dMutRate) + " > " + fnaMutNucs],
                            'targets': [fnaMutNucs],
                            'file_dep': [faaInput],}
                            fnaStartNucs = fnaMutNucs
                        elif strMG == "Illumina-mutated-5pct":
                            fnaMutNucs = dirOutMG + os.sep +strDB + "mutated5pct.fna"
                            yield {'name': "Screen input prot file for nucs, make mutated nucs: " +strDB + strRun + strMG,
                            'actions': ["python " + pyCheckSeqs + " < " + faaInput + " | " + pyMutateNucs + " --mutrate " + ".05" + " > " + fnaMutNucs],
                            'targets': [fnaMutNucs],
                            'file_dep': [faaInput],}
                            fnaStartNucs = fnaMutNucs
                        elif strDB == "ARDB" and strRun!="05" and strMG not in ["Illumina-mutated","Illumina-mutated-5pct"]:
                            fnaStartNucs  = dirData + os.sep + "input_sequences" + os.sep + "ARDBbtnucs.fna"
                        else:
                            fnaStartNucs = dirData + os.sep + "input_sequences" + os.sep + strDB + "nucs.fna"


                        txtScreenLog = dirOutMG + os.sep +strDB + "screeninglog.txt"
                        fnaMGnucs = dirOutMG + os.sep +strDB + "inputMGnucs.fna"


                        yield {'name': "Screen out nucs without a corresponding protein sequence: " +strMG + strRun + strDB,
                       'actions': ["python " + pyScreenNucs + " --out " + fnaMGnucs + " --nucs " + fnaStartNucs + " --prots " + faaCleanProteinSeqs+ " --log " + txtScreenLog ],
                       'targets': [fnaMGnucs],
                       'file_dep': [faaCleanProteinSeqs,fnaStartNucs],}

                        tabProtLengths = dirOutMG + os.sep +strDB + "protlengths.tab"
                        tabNucLengths = dirOutMG + os.sep +strDB + "nuclengths.tab"
                        yield {'name': "Record nuc and prot length: " +tabProtLengths + " " + tabNucLengths,
                       'actions': ["python " + pyCheckSeqs + " < " + faaInput + " |  python " + pyCheckLengths + " > " + tabProtLengths, "python " +  pyCheckLengths + " < " + fnaMGnucs + " > " + tabNucLengths],
                       'targets': [tabProtLengths,tabNucLengths],
                       'file_dep': [faaInput,fnaMGnucs],}

                        dirOutMG = dirSynthMG + os.sep + strDB + os.sep + strMG + "-" + strRun
                        dirFinalGenes = dirOutMG + os.sep + "final_genes"
                        python_check_create_dir(dirSynthMG  + os.sep + strDB)
                        python_check_create_dir(dirSynthMG  + os.sep + strDB + os.sep + strMG)
                        python_check_create_dir(dirOutMG)
                        fastqSim = dirOutMG +os.sep + strMG + "-" + strRun+ "_single.fastq"
                        fastaSim = dirOutMG +os.sep + strMG + "-" + strRun+ ".fasta"

                        # We create the clean genomes for synthetic metagenome separately because GemSim requires all of the material to be in a single folder. A better
                        # approach might be to simply cp these after they finish, or use a softlink.
                        for strGenome in astrGenomes:
                            fnaInGenome = dirBacteriaGenomes + os.sep + strGenome
                            fnaCleanGenome =  dirFinalGenes + os.sep + strGenome
                            dirTmpRemove = dirTmp + os.sep + strDB + "_"+strMG +"_"+ strRun + strGenome
                            python_check_create_dir(dirTmpRemove)
                            tmpGenome = dirTmp + os.sep + strDB + strGenome + strDB + "_"+strMG +"_"+ strRun

                            yield{ 'name': "Remove matching material from genomes: " +strDB + "-" + fnaInGenome + dirOutMG,
                               'actions': ["python " + pyRemovePosGenes + " --prots " +  faaInput + " --genome  " + fnaInGenome + " --tmpgenome " + tmpGenome + " --clean " + fnaCleanGenome + " --tmp " +dirTmpRemove ," rm " + tmpGenome, " rm " + dirTmpRemove + os.sep + "*.*"," rmdir " + dirTmpRemove],
                               'targets': [fnaCleanGenome],
                               'file_dep': [fnaInGenome],}
                            fnaPadGenome = fnaCleanGenome


                        txtGS = dirOutMG + os.sep + strDB + "_"+strMG +"_"+ strRun +"_goldstandard.txt"
                        txtAbundance = dirOutMG + os.sep + strDB + "_"+strMG +"_"+ strRun +"_abundance.txt"
                        txtLogSim = dirOutMG + os.sep + strDB + "_"+strMG +"_"+ strRun +"_log.txt"
                        txtCount = dirOutMG + os.sep + strDB + "_"+strMG +"_"+ strRun +"_readcount.txt"



                        yield{ 'name': "Make Abundance Table and Gold Standard, add material to ends of nucleotide sequences : "+ strDB+"-" +strMG+"-"+strRun,
                                               'actions': ["python " + pySimpleSim + " --nucs " + fnaMGnucs + " -N " + str(iGenes) + " --muS " + str(1) + " --muG " + str(1) + " --gold " + txtGS + " --genomes " + txtGenomes + " --fastadir " + dirFinalGenes + " --abund " + txtAbundance + " --padgenome " + fnaPadGenome + " --padlength "  + str(iPadLength) + " --dirgenomes " +  dirCleanGenomes + os.sep + strDB + " --pctspike " + str(dSpike) + " --log " + txtLogSim],
                                               'targets': [txtAbundance,txtGS],
                                               'file_dep': [fnaMGnucs,fnaPadGenome],}


                        yield{ 'name': "Make synthetic metagenome using GemSim,Convert MG to fasta, : "+ strDB+"-" +strMG+"-"+strRun,
                             'actions': ["python " + pyGemSim + " -R " + dirFinalGenes + " -n " + str(iReads) + " -l d " + " -m " + zipModel + " -c -q " +str(64) + " -o " + dirOutMG + os.sep + strMG + "-" + strRun + " -a " + txtAbundance,"python " + pyFastq2Fasta + " < " + fastqSim + " > " + fastaSim, "rm " + fastqSim ],
                             'targets': [fastaSim],
                             'file_dep': [txtAbundance],}

                        yield{ 'name': "remove excess text in headers, count reads : "+ strDB+"-" +strMG+"-"+strRun,
                             'actions': [" sed -i -e \"s/^>.*\(USR_.*_END\).*$/>\\1/g\" " + fastaSim, "grep -e \">\" -c " + fastaSim + " > " + txtCount],
                             'targets': [txtCount],
                             'file_dep': [fastaSim],}



"""
#############################################################################
Run ShortBRED and centroids on synthetic metagenomes
#############################################################################
* For each clust_id, create ShortBRED markers
(faaInitialMarkers). This will also create corresponding centroids when
ShortBRED clusters the input.
** For each marker length, create a set of ShortBRED markers based on the
data from the initial markers. (This gets around running blastp again.)

* Call SBQUantify for centroids on metagenome.
* Call SBQuantify for markers on metagenome.
"""



adClustID	= [0.80,0.85,0.90,0.95,1.0]
aiAllMarkerLengths = [8,10,12,15,18,20,22,25,30]
c_strRunClust = "85"
c_aiMarkerRun = [8]

c_iMarkerLen = 8
c_dClustID = .85

dirInitialMarkers = dirOut + os.sep + "init_markers"
python_check_create_dir(dirInitialMarkers)


strSBIdentifyLong =  " --qthresh "  + str(c_iQuasiThresh)  + " --xlimit " + str(c_iXLimit) + " --totlength "  + str(c_iTotalMarkerLength) +  " --qclustid " +  str(c_dQuasiClust) + " --len " + str(c_dMatchMaxLen) + " --consthresh " + str(c_dConsThresh) + " --id " + str(c_dMinIDForMatch) + " --qmlength " + str(int(math.floor((100/3)))) +  " --blastp " + c_cmdBlastP + " --muscle " + c_cmdMuscle  + " --cdhit " + c_cmdCDhit + " --usearch " + c_cmdUsearch + " --makeblastdb " + c_cmdMakeBlast


dirSynthMGMarkers = dirOut + os.sep + "synthmg_markers"
dirEvalResults = dirOut + os.sep + "synthmg_results"
tabAllResults = dirEvalResults+  os.sep + "AllSyntheticMGResults.tab"
atabAllStats = []
strStatsHeader = "\t".join(["Dataset","SyntheticMG", "ClustID","MarkerLen","MatchID","SB_AUC","Cent_AUC","SB_Spearman","Cent_Spearman",	"SB_Spec","Cent_Spec","SB_Sens","Cent_Sens"])


def task_EvalSyntheticMetagenomes():
        for strDB in astrDB:
            # Decide how to best replace fnaMutNucs here.
            python_check_create_dir(dirInitialMarkers + os.sep + strDB)

            for dClustID in adClustID:
                strClustID = str(int(dClustID*100))

                # Input for initial marker creation
                faaInput = dirData + os.sep + "input_sequences" + os.sep + strDB + ".input.faa"
                faaClean = dirOut + os.sep + strDB + "clean.faa"

                # Output from initial marker creation
                dirInitMarkers = dirInitialMarkers + os.sep + strDB + strClustID
                faaInitialMarkers = dirInitMarkers + os.sep + strDB + "_" + strClustID + "_initial_Markers.faa"
                faaCentroids = dirInitMarkers + os.sep + "clust" + os.sep + "clust.faa"
                mapCentroids = dirInitMarkers + os.sep + "clust" + os.sep + "clust.map"
                txtBlastRef			= dirInitMarkers + os.sep + "blastresults" + os.sep + "refblast.txt"
                txtBlastSelf			= dirInitMarkers + os.sep + "blastresults" + os.sep + "selfblast.txt"


                yield{ 'name': "Make markers and centroids with ShortBRED-Identify : "+ strDB+"-" + strClustID,
                             'actions': ["python " + pyCheckSeqs + " < " + faaInput +" > " + faaClean , "python " + pySBIdentify  + " --markers " + faaInitialMarkers + " --tmpdir " + dirInitMarkers + " --goi " + faaClean + " --refdb " + c_blastdbRefName + " --threads " + str(c_iBlastThreads) + " --clustid " + str(dClustID) + " --len " + str(c_dMatchMaxLen) + strSBIdentifyParams ],
                             'targets': [faaInitialMarkers],
                             'file_dep': [faaInput],}


                # Prepare folders for markers with different min marker length.
                python_check_create_dir(dirSynthMGMarkers)
                python_check_create_dir(dirSynthMGMarkers + os.sep + strDB + "_"+ strClustID)
                python_check_create_dir(dirEvalResults)
                python_check_create_dir(dirEvalResults+  os.sep + strDB )

                for iMarkerLen in aiAllMarkerLengths:
                    dirTmpMarkers = dirSynthMGMarkers + os.sep + strDB +"_"+ strClustID + os.sep +  "MarkerTmp"+str(iMarkerLen)
                    python_check_create_dir(dirTmpMarkers)
                    faaMarkers = dirTmpMarkers + os.sep + strDB+str(iMarkerLen) + "markers.faa"

                    yield{ 'name': "Make additional markers for each marker length : "+ strDB+"-" + strClustID+ "-" + str(iMarkerLen),
                             'actions': ["python " + pySBIdentify  + " --markers " + faaMarkers + " --tmpdir " + dirTmpMarkers + " --goiclust " + faaCentroids + " --map_in " + mapCentroids + " --refblast " + txtBlastRef + " --goiblast " + txtBlastSelf + " --markerlength " + str(iMarkerLen) + strSBIdentifyLong ],
                             'targets': [faaMarkers],
                             'file_dep': [faaInitialMarkers,faaCentroids,mapCentroids,txtBlastRef,txtBlastSelf]}

                # Run centroids and markers through ShortBRED Quantify
                for strMG in ["Illumina","Illumina-mutated","Illumina-mutated-5pct"]:
                    for strRun in ["05","10","25"]:
                        dirCentEval = dirEvalResults+  os.sep + strDB + os.sep + strMG +"_"+ strRun + os.sep + strClustID + os.sep + "cent_results"
                        python_check_create_dir(dirEvalResults+  os.sep + strDB + os.sep + strMG +"_"+ strRun)
                        python_check_create_dir(dirEvalResults+  os.sep + strDB + os.sep + strMG +"_"+ strRun + os.sep + strClustID)
                        python_check_create_dir(dirEvalResults+  os.sep + strDB + os.sep + strMG +"_"+ strRun + os.sep + strClustID + os.sep + "cent_results")

                        dirOutMG = dirSynthMG + os.sep + strDB + os.sep + strMG + "-" + strRun
                        dirMGData = dirSynthMG + os.sep + strDB + os.sep + strMG + "-" + strRun
                        fnaSimData = dirMGData +os.sep + strMG + "-" + strRun+ ".fasta"



                        fnaMGnucs = dirOutMG + os.sep +strDB + "inputMGnucs.fna"

                        txtGS = dirOutMG + os.sep + strDB + "_"+strMG +"_"+ strRun +"_goldstandard.txt"
                        txtAbundance = dirOutMG + os.sep + strDB + "_"+strMG +"_"+ strRun +"_abundance.txt"
                        txtLogSim = dirOutMG + os.sep + strDB + "_"+strMG +"_"+ strRun +"_log.txt"
                        txtCount = dirOutMG + os.sep + strDB + "_"+strMG +"_"+ strRun +"_readcount.txt"
                        tabProtLengths = dirOutMG + os.sep +strDB + "protlengths.tab"

                        blastCentroids = dirCentEval + os.sep + "blastCentroids.tab"
                        blastFullHits =  dirCentEval + os.sep + "blastCentroids_fullhits.tab"
                        txtCentResults = dirCentEval + os.sep + "centroid_results.tab"
                        txtTimeCent =  dirCentEval + os.sep + "time.txt"

                        # Run all possible clustering levels and marker lengths on the 5% Illumina metagenome.
                        if strRun=="05" and strMG=="Illumina":
                            aiMarkerLenMGRun = aiAllMarkerLengths
                            bRunMarkersAndCentroids = True

                        # On other metagenomes, only run one choice of marker length and clustering identity.
                        # The marker length and clustering id were chosen based on the results for the 5pct metagenome.
                        # This was split more distinctively in the original shortbred-sfle project, and is combined here
                        # to make rerunning easier.
                        else:
                            aiMarkerLenMGRun = c_aiMarkerRun
                            bRunMarkersAndCentroids = (strClustID == c_strRunClust)

                        if(bRunMarkersAndCentroids):
                            yield{ 'name': "Run centroids on metagenome:" + " ".join([strDB,strClustID,strMG,strRun]),
                            'actions':["/usr/bin/time "+" -o " + txtTimeCent +" "+pySBQuantify+" --markers "+faaCentroids+" --threads 5 --wgs "+fnaSimData+" --tmp "+dirCentEval+" --blastout "+blastCentroids+" --notmarkers "+" Y "+" --id "+ str(c_dMatchID)+" --SBhits "+blastFullHits+" --results "+txtCentResults+" --small "+"True" + " --usearch "+ c_cmdUsearch],
                            'targets':[txtCentResults,blastCentroids,blastFullHits], 'file_dep':[faaCentroids,fnaSimData]}


                            for iMarkerLen in aiMarkerLenMGRun:
                                dirTmpMarkers = dirSynthMGMarkers + os.sep + strDB +"_"+ strClustID + os.sep +  "MarkerTmp"+str(iMarkerLen)
                                faaMarkers = dirTmpMarkers + os.sep + strDB+str(iMarkerLen) + "markers.faa"
                                dirMarkerEval = dirEvalResults+  os.sep + strDB + os.sep + strMG +"_"+ strRun + os.sep + strClustID + os.sep + "markers" + str(iMarkerLen)
                                python_check_create_dir(dirMarkerEval)

                                blastMarkers = dirMarkerEval + os.sep + "blastMarkers.tab"
                                blastMarkerFullHits =  dirMarkerEval + os.sep + "blastMarkers_fullhits.tab"
                                txtMarkerResults = dirMarkerEval + os.sep + "marker_results.tab"
                                txtMarkerOut = dirMarkerEval + os.sep + "markersUpdatedRPKM.tab"
                                txtTimeMarker =  dirMarkerEval + os.sep + "time.txt"

                                csvWGSCounts = dirMarkerEval + os.sep + "gene_counts.csv"
                                mapFinal = dirTmpMarkers+ os.sep + "final.map"
                                iReadCount = 5000000
                                txtMarkerFP = dirMarkerEval + os.sep + "MarkerFPs.tab"
                                tabCGResults = dirMarkerEval + os.sep + "prelim_eval_results.tab"
                                tabLow = dirMarkerEval + os.sep + "low.tab"
                                tabHigh = dirMarkerEval + os.sep + "high.tab"
                                txtEvalTable =  dirMarkerEval + os.sep + "full_eval_results.tab"


                                yield{ 'name': "Run markers on metagenome:" + " ".join([strDB,strClustID,strMG,strRun,str(iMarkerLen)]),
                            'actions':["/usr/bin/time "+" -o " + txtTimeMarker +" "+pySBQuantify+" --markers "+faaMarkers+" --threads 5 --wgs "+fnaSimData+" --tmp "+dirMarkerEval+" --blastout "+blastMarkers+" --id "+ str(c_dMatchID)+" --SBhits "+blastMarkerFullHits+" --results "+txtMarkerResults+" --marker_results " + txtMarkerOut +" --small "+"True" + " --usearch "+ c_cmdUsearch],
                            'targets':[txtMarkerResults,blastMarkers,blastMarkerFullHits,txtMarkerOut], 'file_dep':[faaMarkers,fnaSimData]}

                                yield{ 'name': "Compile centroid and marker results for evaluation:" + " ".join([strDB,strClustID,strMG,strRun,str(iMarkerLen)]),
                            'actions':["python " + pyCountWGS + " --map_in " + mapFinal + " --wgs " + fnaSimData + " > " + csvWGSCounts, "python " + pyCheckGenes + " --reads " + str(iReadCount) +	" --wgs " + csvWGSCounts + " --markers " + faaMarkers + " --centroids " + faaCentroids + " --markerblast " + blastMarkerFullHits + " --centblast " + blastFullHits + " --qmarkers " + txtMarkerResults + " --qcent "  + txtCentResults + " --map " + mapFinal + " --nucs " + fnaMGnucs + " --FPfile " + txtMarkerFP + " > " + tabCGResults],
                                       'targets':[tabCGResults,txtMarkerFP,csvWGSCounts],'file_dep':[txtMarkerOut,fnaMGnucs,txtCentResults]}



                                yield{ 'name': "Merge on additional evaluation files:" + " ".join([strDB,strClustID,strMG,strRun,str(iMarkerLen)]),
                            'actions':["Rscript "  + rMergeData + " --gs " + txtGS +  " --wc " + csvWGSCounts + " --main " + tabCGResults + " --out " + txtEvalTable + " --cent " + tabProtLengths + " --low " + tabLow + " --high " +  tabHigh + " --markers " + txtMarkerOut],
                                       'targets':[txtEvalTable],'file_dep':[csvWGSCounts,tabCGResults]}

                                pngSpec	   		= dirMarkerEval + os.sep + "spec.png"
                                pngCorrelation	  	= dirMarkerEval + os.sep + "corr.png"
                                pngTprFpr	  		= dirMarkerEval + os.sep + "tprXfpr.png"
                                txtStats	   		= dirMarkerEval + os.sep + strDB+"p"+strClustID + "m"+str(iMarkerLen)+".tab"
                                atabAllStats.append(txtStats)

                                yield{ 'name': "Plot each run, generate summary stats:" + " ".join([strDB,strClustID,strMG,strRun,str(iMarkerLen)]),
                            'actions':["Rscript "  + rEval + " --name " + strDB + " --rf " + txtEvalTable + " --corr " + pngCorrelation + " --spec " + pngSpec + " --tpfp " + pngTprFpr + " --stats " + txtStats + " --ml " + str(iMarkerLen) + " --pct " +  strClustID + " --id " + str(c_dMatchID) + " --metagenome " + strMG+"_"+strRun],
                                       'targets':[pngSpec,pngCorrelation,pngTprFpr,txtStats],'file_dep':[txtEvalTable]}

        yield{ 'name': "Combine all statistics:" + " ".join([strDB,strClustID,strMG,strRun]),
    'actions':["echo \" " + strStatsHeader + " \" " + " | cat  - " + " ".join(atabAllStats) + " > " + tabAllResults,"Rscript " + rPlotAllAUC + " --AUC " + tabAllResults + " --db VF --outdir " + dirFigures, "Rscript " + rPlotAllAUC + " --AUC " + tabAllResults + " --db ARDB --outdir " + dirFigures ],
           'targets':[tabAllResults],'file_dep':atabAllStats}


        astrMGFig2Results = ["/n/huttenhower_lab/data/shortbred/shortbred_doit/out/synthmg_results/ARDB/Illumina_05/85/markers8/full_eval_results.tab",
                              "/n/huttenhower_lab/data/shortbred/shortbred_doit/out/synthmg_results/ARDB/Illumina_10/85/markers8/full_eval_results.tab",
                              "/n/huttenhower_lab/data/shortbred/shortbred_doit/out/synthmg_results/ARDB/Illumina_25/85/markers8/full_eval_results.tab",
                              "/n/huttenhower_lab/data/shortbred/shortbred_doit/out/synthmg_results/VF/Illumina_05/85/markers8/full_eval_results.tab",
                              "/n/huttenhower_lab/data/shortbred/shortbred_doit/out/synthmg_results/VF/Illumina_10/85/markers8/full_eval_results.tab",
                              "/n/huttenhower_lab/data/shortbred/shortbred_doit/out/synthmg_results/VF/Illumina_25/85/markers8/full_eval_results.tab"]

        # Produce Figure 2
        pdfFig2 = dirFigures + os.sep + "Fig2.pdf"
        epsFig2 = dirFigures + os.sep + "Fig2.eps"
        yield{ 'name': "Figure 2:",
        'actions':["Rscript " + rFig2 + " --A05  " + astrMGFig2Results[0] + " --A10  " + astrMGFig2Results[1] + " --A25  " + astrMGFig2Results[2] +
      " --V05  " + astrMGFig2Results[3] + " --V10  " + astrMGFig2Results[4] + " --V25  " + astrMGFig2Results[5] + " --fig " + pdfFig2 + " --AUC " +
      dirFigures +os.sep + "Fig2.tab" + " --reads 5000000 --title \"\" --fullROC Y --eps " + epsFig2],
      'targets':[epsFig2],'file_dep':astrMGFig2Results,'uptodate':[False]}


#	oo.pipe([atxtEvalResults],[pdfFig2,tabAUC_MG,epsFig2],RFig2, A05=strA05,
#	A10=strA10, A25=strA25, V05=strV05, V10=strV10, V25=strV25,fig=pdfFig2 ,eps=epsFig2,reads=iReadCount,title="",fullROC=strFullROC,AUC=tabAUC_MG)
#	Default([pdfFig2,tabAUC_MG,epsFig2])


"""
#############################################################################
Build Markers and Centroids for yeast, test against a bacterial metagenome
#############################################################################



fnaMetagenome =  ditOut + os.sep + "synthetic_metagenomes/ARDB/Illumina-05/Illumina-05.fasta"

# Constants
iLimitProts = 2000
dOrgClustID = 0.85
iOrgThreads = 100
iMarkerLen = 8
strOrgDir = "testorg"
iOrgQuasiThresh = c_iQuasiThresh
iOrgMarkerLength = 8


pyLimitSeqs = dirSrc + os.sep +"utils"+ os.sep + "LimitSeqs.py"

faaInputProts = dirData+os.sep+"yeast"+"*.pep"
strOrganism = os.path.basename(str(faaInputProts)).replace( ".pep", "" )
print strOrganism

#--Tmp files for Cleaning Org Prots--#
dirYeast =dirTmp+os.sep +"yeast"
python_check_create_dir(dirYeast)
faaOrgProts = faaInputProts
fnaOrgNucs = faaInputProts.replace(".pep",".genome")
faaOrgCleanProts = dirYeast + os.sep + strOrganism+"-clean.pep"
faaOrgReducedProts = dirYeast + os.sep + strOrganism+"-small.pep"

dirOrgMarkers = diryYeast + os.sep + "SBIdentify"

faaOrgInitialMarkers = dirOrgMarkers + os.sep + "yeast_markers.faa"
faaOrgCentroids  =  dirOrgMarkers + os.sep + "clust" + "clust.faa"

yield{ 'name': "BuildMarkersAndCentroids for Yeast:" + faaInputProts),
'actions':["python " + pyCheckSeqs + " < " + faaInputProts +" > " + faaOrgCleanProts , "python " + pyLimitSeqs + " -N " + str(iLimitProts) + " < " + faaOrgCleanProts + " > " + faaOrgReducedProts], "python " + pySBIdentify  + " --markers " + faaOrgInitialMarkers + " --tmpdir " + dirOrgMarkers + " --goi " + faaOrgReducedProts + " --refdb " + c_blastdbRefName + " --threads " + str(c_iBlastThreads) + " --clustid " + str(c_dClustID) + " --len " + str(c_dMatchMaxLen) + strSBIdentifyParams + " --markerlength " + str(c_iMarkerLen) + " --qthresh " + str(c_iQuasiThresh)],'targets':[faaOrgInitialMarkers],'file_dep':[faaInputProts]}


    #Run Shortbred Markers against the wgs data
   	oo.pipe([faaOrgMarkers,fnaMetagenome],[txtOrgMarkerResults],ShortbredQ,
	markers=faaOrgMarkers,threads=iQuantThreads,wgs=fnaMetagenome,tmp=dirOrgMarkerQuant,blastout=blastOrgMarkers,
	id=c_dMatchID,
	results=txtOrgMarkerResults)
	Default(txtOrgMarkerResults)



    #Run centroids against the wgs data
	oo.pipe([faaOrgCentroids,fnaMetagenome],[txtOrgCentResults],ShortbredQ,
	markers=faaOrgCentroids,threads=iQuantThreads,wgs=fnaMetagenome,tmp=dirOrgCentQuant,blastout=blastOrgCentroids,
	id=c_dMatchID,notmarkers="Y",results=txtOrgCentResults)

	Default(txtOrgCentResults,txtOrgMarkerResults)
"""

"""
#############################################################################
Build ShortBRED Markers
#############################################################################
* Create Initial Markers for ARDBfilter_WashU.input.faa
* Use the data from initial run to create markers for HMP, T2D_Short, and T2D_long
and Bacteria
"""

dARClustID = .85
strARClustID = str(int(dARClustID*100))
strAR_DB = "ARDBfilter_WashU"
# Input for initial marker creation
faaARInput = dirData + os.sep + "input_sequences" + os.sep + "ARDBfilter_WashU.input.faa"
faaARClean = dirOut + os.sep + "init_markers" + os.sep + strAR_DB + "clean.faa"
dirARInitMarkers = dirOut + os.sep + "init_markers" + os.sep + strAR_DB + os.sep + strARClustID
# Output from initial marker creation
faaARInitialMarkers = dirARInitMarkers + os.sep + strAR_DB + "_" + strARClustID + "_initial_Markers.faa"
faaARCentroids = dirARInitMarkers + os.sep + "clust" + os.sep + "clust.faa"
mapARCentroids = dirARInitMarkers + os.sep + "clust" + os.sep + "clust.map"
txtARBlastRef			= dirARInitMarkers + os.sep + "blastresults" + os.sep + "refblast.txt"
txtARBlastSelf			= dirARInitMarkers + os.sep + "blastresults" + os.sep + "selfblast.txt"
python_check_create_dir(dirOut + os.sep + "init_markers" + os.sep + strAR_DB)
python_check_create_dir(dirARInitMarkers)

def task_MakeInitial_AR_Markers():
    return {
        'actions': ["python " + pyCheckSeqs + " < " + faaARInput +" > " + faaARClean , "python " + pySBIdentify  + " --markers " + faaARInitialMarkers + " --tmpdir " + dirARInitMarkers + " --goi " + faaARClean + " --refdb " + c_blastdbRefName + " --threads " + str(150) + " --clustid " + str(dARClustID) + " --len " + str(c_dMatchMaxLen) + strSBIdentifyParams ],
        'file_dep': [faaARInput],
        'targets': [faaARInitialMarkers],
    }

dARClustID = .85
strARClustID = str(int(dARClustID*100))
strAR_DB = "ARDBfilter_WashU"
# Input for initial marker creation
faaARInput = dirData + os.sep + "input_sequences" + os.sep + "ARDBfilter_WashU.input.faa"
faaARClean = dirOut + os.sep + "init_markers" + os.sep + strAR_DB + "clean.faa"

dirARSandyInitMarkers = dirOut + os.sep + "init_markers_sandy" + os.sep + strAR_DB + os.sep + strARClustID
# Output from initial marker creation
faaARSandyInitialMarkers = dirARSandyInitMarkers + os.sep + strAR_DB + "_" + strARClustID + "sandy_initial_Markers.faa"
faaARSandyCentroids = dirARSandyInitMarkers + os.sep + "clust" + os.sep + "clust.faa"
mapARSandyCentroids = dirARSandyInitMarkers + os.sep + "clust" + os.sep + "clust.map"
txtARSandyBlastRef			= dirARSandyInitMarkers + os.sep + "blastresults" + os.sep + "refblast.txt"
txtARSandyBlastSelf			= dirARSandyInitMarkers + os.sep + "blastresults" + os.sep + "selfblast.txt"
python_check_create_dir(dirOut + os.sep + "init_markers_sandy")
python_check_create_dir(dirOut + os.sep + "init_markers_sandy" + os.sep + strAR_DB)
python_check_create_dir(dirARSandyInitMarkers)


def task_MakeInitial_AR_Markers_SandyParams():
    return {
        'actions': ["python " + pySBIdentify  + " --markers " + faaARSandyInitialMarkers + " --tmpdir " +
        dirARSandyInitMarkers + " --goi " + faaARClean + " --refdb " + c_blastdbRefName + " --threads " + str(150) +
        " --clustid " + str(dARClustID) + " --len " + str(c_dMatchMaxLen) +  " --consthresh " +str(c_dConsThresh) +
        " --blastp " + c_cmdSandyBlastP + " --muscle " + c_cmdSandyMuscle  + " --cdhit " + c_cmdSandyCDhit + " --usearch " + c_cmdUsearch +
        " --makeblastdb " + c_cmdSandyMakeBlast ],
        'file_dep': [faaARClean],
        'targets': [faaARSandyInitialMarkers],
    }
#oo.pipe([faaCleanWashUProts, c_blastdbRefPal], [faaInitialWashUMarkers,txtWashUBlastRef,
#	txtWashUBlastSelf,mapWashUCentroids,faaWashUCentroids,mapWashU], ShortbredID,
#	markers=faaInitialWashUMarkers,
#	tmpdir=tmpDirectory ,goi=faaCleanWashUProts, refdb=c_blastdbRefName,
#	threads=150, clustid=dPercentIdentity,len=c_dMatchMaxLen,
#	consthresh=c_dConsThresh,blastp=c_cmdBlastP)

astrMarkerSets = ["HMP","T2D","T2D_short","Bacteria"]
dirHumanMGMarkers = dirOut + os.sep + "human_mg_markers"
python_check_create_dir(dirHumanMGMarkers)


def task_MakeAdditionalARMarkers():
    for strMG in astrMarkerSets:
        if strMG == "HMP":
            iAvgBP = 101
            iMinTrustBP = 90
        elif strMG == "T2D":
            iAvgBP = 90
            iMinTrustBP = 81
        elif strMG == "T2D_short":
            iAvgBP = 75
            iMinTrustBP = 68
        elif strMG == "Bacteria":
            iAvgBP = 100
            iMinTrustBP = 30

        faaMGMarkers = dirHumanMGMarkers + os.sep + strMG + ".faa"
        faaSandyMGMarkers = dirHumanMGMarkers + os.sep + strMG + "_sandy.faa"
        dirTmpMG = dirHumanMGMarkers + os.sep + "SBdata_" + strMG
        dirTmpSandyMG = dirHumanMGMarkers + os.sep + "SBdata_Sandy" + strMG


        yield{ 'name': "Make additional markers for each set : "+ strMG,
                             'actions': ["python " + pySBIdentify  + " --markers " + faaMGMarkers + " --tmpdir " + dirTmpMG + " --goiclust " + faaARCentroids + " --map_in " + mapARCentroids + " --refblast " + txtARBlastRef + " --goiblast " + txtARBlastSelf + " --markerlength " + str(c_iMarkerLen) + " --qthresh " + str(c_iQuasiThresh) + " --qmlength " + str(int(math.floor((iAvgBP/3)))) + " --qclustid " + str(c_dQuasiClust) + " --len " + str(c_dMatchMaxLen) + " --consthresh " + str(c_dConsThresh) + " --id " + str(c_dMinIDForMatch) + " --totlength " + str(300) + " --muscle " + c_cmdMuscle  + " --cdhit " + c_cmdCDhit + " --usearch " + c_cmdUsearch + " --blastp " + c_cmdBlastP + " --makeblastdb " + c_cmdMakeBlast ],
                             'targets': [faaMGMarkers],
                             'file_dep': [faaARInitialMarkers,faaARCentroids,mapARCentroids,txtARBlastRef,txtARBlastSelf]}
        yield{ 'name': "Make additional sandy markers for each set : "+ faaSandyMGMarkers,
                             'actions': ["python " + pySBIdentify  + " --markers " + faaSandyMGMarkers + " --tmpdir " + dirTmpSandyMG + " --goiclust " + faaARSandyCentroids +
                             " --map_in " + mapARSandyCentroids + " --refblast " + txtARSandyBlastRef + " --goiblast " + txtARSandyBlastSelf + " --markerlength " + str(c_iMarkerLen) +
                             " --qthresh " + str(c_iQuasiThresh) + " --qmlength " + str(int(math.floor((iAvgBP/3)))) + " --qclustid " + str(c_dQuasiClust) +
                             " --len " + str(c_dMatchMaxLen) + " --consthresh " + str(c_dConsThresh) + " --id " + str(c_dMinIDForMatch) + " --totlength " + str(300) +
                             " --muscle " + c_cmdSandyMuscle  + " --cdhit " + c_cmdSandyCDhit + " --usearch " + c_cmdUsearch + " --blastp " + c_cmdSandyBlastP + " --makeblastdb " + c_cmdSandyMakeBlast ],
                             'targets': [faaSandyMGMarkers],
                             'file_dep': [faaARSandyInitialMarkers,faaARSandyCentroids,mapARSandyCentroids,txtARSandyBlastRef,txtARSandyBlastSelf]}

#oo.pipe([txtBlastRef,txtBlastSelf,mapCentroids,faaCentroids], [faaWGSMarkers, mapWGSFinal], ShortbredID,
#					markers=faaWGSMarkers, tmpdir=dirMarker ,goiclust=faaCentroids, map_in=mapCentroids,
#					refblast=txtBlastRef, goiblast=txtBlastSelf, threads=iSBThreads+1, clustid=dPercentIdentity,
#					qthresh=c_iQuasiThresh,markerlength=iMarkerLen,qmlength=int(math.floor((iAvgBP/3))),
#					qclustid=c_dQuasiClust,len=c_dMatchMaxLen,consthresh=c_dConsThresh,id=c_dMinIDForMatch,totlength=300)


"""
#############################################################################
Profile Human Metagenomes with ShortBRED
#############################################################################
* Compile the list of samples in dataset.
* Profile metagenomes in dataset.
* Merge Results
"""

"""
dirMGInfo = dirData +os.sep + "human_mg_info"
txtHMPSampleList = dirMGInfo + os.sep + "TP1_and_QC.txt"
txtT2DSampleList = dirMGInfo + os.sep + "T2D_list.txt"

dirHMPSamples = "/n/huttenhower_lab_nobackup/downloads/HMP/TP1_stool/"
dirT2DSamples = "/n/huttenhower_lab_nobackup/downloads/T2D/samples/fastq/"
dirYATSamples = "/n/huttenhower_lab_nobackup/downloads/Yatsunenko/wgs/"

#Load the T2D Read Lengths into a dictionary. They vary from sample to sample.
txtT2DReadLengths = dirMGInfo + os.sep  + "T2D_readlengths.txt"
dictT2DReadLengths = {}

dirMGResults = dirOut + os.sep + "human_mg_results"
python_check_create_dir(dirMGResults)

def task_RunShortBREDOnHumanMetagenomes():

    for strWGSdb in ["YAT","T2D","T2D_short","HMP"]:
        astrWGSSamples = []
        astrResultsFiles = []
        dirWGSResults 	= dirMGResults + os.sep + strWGSdb
        python_check_create_dir(dirWGSResults)
        txtMerged = dirWGSResults + os.sep +strWGSdb+"mergedresults.tab"

        if strWGSdb == "YAT":
            dirWGSSamples = dirYATSamples
            afileWGSSamples = glob.glob(dirWGSSamples + os.sep + "*.fna")
            YATResults = txtMerged
            for fileSample in afileWGSSamples:
                astrWGSSamples.append(str(fileSample))
        elif strWGSdb == "HMP":
            dirWGSSamples = dirHMPSamples
            HMPResults = txtMerged
            for strLine in open(txtHMPSampleList,'r'):
                strLine = dirWGSSamples + strLine.strip()
                #strLine = dirWGSSamples + strLine.strip() + ".tar.bz2"
                astrWGSSamples.append(strLine)
        elif strWGSdb == "T2D":
            dirWGSSamples = dirT2DSamples
            T2DResults = txtMerged
            for strLine in open(txtT2DSampleList,'r'):
                astrLine = strLine.split(" ")
                strLine1 = astrLine[0].strip()
                strLine2 = astrLine[1].strip()
                mtchID = re.search(r'_n_CHB_data_T2D_samples_(.*)\.sra_1\.fastq\.gz',astrLine[0])
                # Not all read lengths are included in the T2D metadata. I assume 75 bp for each undefined one.
                iFileBP = dictT2DReadLengths.get(mtchID.group(1),75)
                print(mtchID.group(1),iFileBP)
                if iFileBP == 90:
                    #print("Added File!",strLine1)
                    astrWGSSamples.append([strLine1,strLine2])
        elif strWGSdb == "T2D_short":
            dirWGSSamples = dirT2DSamples
            T2DResults = txtMerged
            for strLine in open(txtT2DSampleList,'r'):
                astrLine = strLine.split(" ")
                strLine1 = astrLine[0].strip()
                strLine2 = astrLine[1].strip()
                mtchID = re.search(r'_n_CHB_data_T2D_samples_(.*)\.sra_1\.fastq\.gz',astrLine[0])
                # Not all read lengths are included in the T2D metadata. I assume 75 bp for each undefined one.
                iFileBP = dictT2DReadLengths.get(mtchID.group(1),75)
                if iFileBP == 75:
					astrWGSSamples.append([strLine1,strLine2])

        print astrWGSSamples

        for strSample in astrWGSSamples:
            if strWGSdb == "YAT":
                strOutName = dirWGSResults + os.sep + os.path.basename(strSample).replace(".fna","")+"RPKM.out"
                dirQuant = c_dirBigTmp + os.sep + strWGSdb + os.sep + os.path.basename(strSample).replace(".fna","")
                iAvgBP = 450
                iMinTrustBP = 200
            elif strWGSdb == "HMP":
                strOutName = dirWGSResults + os.sep + os.path.basename(strSample).replace(".tar.bz2","")+"RPKM.out"
                dirQuant = c_dirBigTmp + os.sep + strWGSdb + os.sep + os.path.basename(strSample).replace(".tar.bz2","")
                iAvgBP = 101
                iMinTrustBP = 90
            elif strWGSdb == "T2D":
                astrSample = strSample
                mtchID = re.search(r'_n_CHB_data_T2D_samples_(.*)\.sra_1\.fastq\.gz',astrSample[0])
                strSample = dirT2DSamples + astrSample[0] + " " +dirT2DSamples + astrSample[1]
                sys.stderr.write(strSample + "\n")
                strOutName = dirWGSResults + os.sep + mtchID.group(1)+"RPKM.out"
                dirQuant = c_dirBigTmp + os.sep + strWGSdb + os.sep + mtchID.group(1)
                iAvgBP = dictT2DReadLengths.get(mtchID.group(1),90)
                iMinTrustBP = (float(iAvgBP) * .90)
            elif strWGSdb == "T2D_short":
                astrSample = strSample
                mtchID = re.search(r'_n_CHB_data_T2D_samples_(.*)\.sra_1\.fastq\.gz',astrSample[0])
                strSample = dirT2DSamples + astrSample[0] + " " +dirT2DSamples + astrSample[1]
                sys.stderr.write(strSample + "\n")
                strOutName = dirWGSResults + os.sep + mtchID.group(1)+"RPKM.out"
                dirQuant = c_dirBigTmp +  os.sep + strWGSdb + os.sep + mtchID.group(1)
                iAvgBP = dictT2DReadLengths.get(mtchID.group(1),75)
                iMinTrustBP = (float(iAvgBP) * .90)

            txtCleanResults = strOutName.replace(".out",".txt")
            if strWGSdb == "YAT":
                faaMGMarkers = faaARCentroids
                strCent = "Y"
            else:
                faaMGMarkers = dirHumanMGMarkers + os.sep + strWGSdb+ ".faa"
                strCent = "N"

            print "python " + pySBQuantify+" --markers "+faaMGMarkers+" --threads 2 --wgs "+strSample+" --tmp "+dirQuant+" --id "+ str(c_dMatchID)+" --results "+strOutName+" --minreadBP " + str(iMinTrustBP) + " --avgreadBP " + str(iAvgBP) + " --usearch "+ c_cmdUsearch + " --notmarkers " + strCent

            yield{ 'name': "Run markers on human metagenome:" + strWGSdb + "-" + strSample,
            'actions':["python " + pySBQuantify+" --markers "+faaMGMarkers+" --threads 5 --wgs "+strSample+" --tmp "+dirQuant+" --id "+ str(c_dMatchID)+" --results "+strOutName+" --minreadBP " + str(iMinTrustBP) + " --avgreadBP " + str(iAvgBP) + " --usearch "+ c_cmdUsearch + " --notmarkers " + strCent, "cut " + strOutName  + " -f 1,2 > "+ txtCleanResults],'targets':[strOutName], 'file_dep':[faaMGMarkers]}
"""




