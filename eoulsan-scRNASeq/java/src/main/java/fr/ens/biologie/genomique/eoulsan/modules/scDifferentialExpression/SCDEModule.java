package fr.ens.biologie.genomique.eoulsan.modules.scDifferentialExpression;

import fr.ens.biologie.genomique.eoulsan.EoulsanException;
import fr.ens.biologie.genomique.eoulsan.Globals;
import fr.ens.biologie.genomique.eoulsan.annotations.LocalOnly;
import fr.ens.biologie.genomique.eoulsan.core.*;
import fr.ens.biologie.genomique.eoulsan.core.workflow.TaskContextImpl;
import fr.ens.biologie.genomique.eoulsan.data.*;
import fr.ens.biologie.genomique.eoulsan.design.Design;
import fr.ens.biologie.genomique.eoulsan.design.Experiment;
import fr.ens.biologie.genomique.eoulsan.design.ExperimentMetadata;
import fr.ens.biologie.genomique.eoulsan.galaxytools.executorinterpreters.DockerExecutorInterpreter;
import fr.ens.biologie.genomique.eoulsan.galaxytools.executorinterpreters.ExecutorInterpreter;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import static fr.ens.biologie.genomique.eoulsan.EoulsanLogger.getLogger;

import static fr.ens.biologie.genomique.eoulsan.core.OutputPortsBuilder.singleOutputPort;

/**
 * Created by brelurut on 18/01/17.
 */

@LocalOnly public class SCDEModule extends AbstractSCDifferentialExpression {

    /**
     * Module Name
     */
    private static final String MODULE_NAME = "scde";

    /**
     * R executor
     */
    private ExecutorInterpreter executor;

    /**
     * Base command
     */
    private static final String cmd = "Rscript";

    /**
     * Scripts position
     */
    private static final String modelScript = "/scripts/SCEscdeErrors.R";
    private static final String priorScript = "/scripts/SCEscdePriors.R";
    private static final String testScript = "/scripts/SCEscdeTest.R";

    /**
     * Module output prefixes
     */
    private static final String MODULE_PREFIX = "scde";
    private static final String MODEL_PREFIX = MODULE_PREFIX + "ModelMatrix";
    private static final String PRIOR_PREFIX = MODULE_PREFIX + "PriorModel";

    /**
     * Parameters Names
     */

    //
    // Main Parameters
    //

    private static final String DOCKER_IMAGE = "docker";

    private static final String MODEL_FIT_COL = "model.group.col";

    private static final String PRIOR_LENGTH = "prior.length";

    private static final String BATCH_COL = "batch.col";
    private static final String RANDOMIZATIONS = "n.randomizations";

    //
    // More Options
    //

    // Model fitting Options
    private static final String MIN_OBS = "min.observation";
    private static final String MIN_GENES = "min.genes";

    private static final String THRESHOLD_SEG = "threshold.segmentation";
    private static final String FAILURE_THRESHOLD = "failure.threshold";

    private static final String MAX_PAIRS = "max.pairs";
    private static final String MIN_PAIRS = "min.pairs";

    private static final String POISSON_PARAM = "poisson.param";
    private static final String LINEAR_FIT = "linear.fit";
    private static final String MIN_THETA = "min.theta";
    private static final String MAX_THETA = "max.theta";

    private static final String MODEL_PLOTS = "save.model.plots";

    // Prior calculation Options
    private static final String PRIOR_PLOT = "save.prior.plot";
    private static final String PSEUDOCOUNT = "pseudocount";
    private static final String QUANTILE = "quantile";
    private static final String MAX_VALUE = "max.value";

    // Test Options
    private static final String POSTERIORS = "return.posteriors";

    /**
     * Default Parameters
     */

    private static final String DEFAULT_PLOTS = "TRUE";

    private static final String DEFAULT_DOCKER =
            "genomicpariscentre/dge-scde:3.7";

    private String docker = DEFAULT_DOCKER;

    // Model fitting parameters and options

    private static final String DEFAULT_GROUP_COL = "NULL";

    private static final int DEFAULT_MIN_OBS = 3;
    private static final int DEFAULT_MIN_GENES = 2000;

    private static final String DEFAULT_THRESHOLD_SEG = "TRUE";
    private static final int DEFAULT_FAILURE_THRESHOLD = 4;

    private static final int DEFAULT_MAX_PAIRS = 5000;
    private static final int DEFAULT_MIN_PAIRS = 10;

    private static final double DEFAULT_POISSON_PARAM = 0.1;

    private static final String DEFAULT_LINEAR_FIT = "TRUE";
    private static final double DEFAULT_MIN_THETA = 1e-2;
    private static final double DEFAULT_MAX_THETA = 1e2;

    private String groupCol = DEFAULT_GROUP_COL;

    private int minObservation = DEFAULT_MIN_OBS;
    private int minGenes = DEFAULT_MIN_GENES;

    private String thresholdSeg = DEFAULT_THRESHOLD_SEG;
    private int failureThreshold = DEFAULT_FAILURE_THRESHOLD;

    private int maxPairs = DEFAULT_MAX_PAIRS;
    private int minPairs = DEFAULT_MIN_PAIRS;

    private double poissonParam = DEFAULT_POISSON_PARAM;

    private String linearFit = DEFAULT_LINEAR_FIT;
    private double minTheta = DEFAULT_MIN_THETA;
    private double maxTheta = DEFAULT_MAX_THETA;

    private String modelPlot = DEFAULT_PLOTS;

    // Priors calculation parameters and options

    private static final int DEFAULT_PRIOR_LENGTH = 400;
    private static final int DEFAULT_PSEUDOCOUNT = 1;
    private static final double DEFAULT_QUANTILE = 0.999;
    private static final Integer DEFAULT_MAX_VALUE = null;

    private int length = DEFAULT_PRIOR_LENGTH;
    private String priorPlot = DEFAULT_PLOTS;
    private int pseudocount = DEFAULT_PSEUDOCOUNT;
    private double quantile = DEFAULT_QUANTILE;
    private Integer maxValue = DEFAULT_MAX_VALUE;

    // Differential expression parameters and options

    private static final int DEFAULT_RANDOMIZATIONS = 150;
    private static final String DEFAULT_POSTERIORS = "TRUE";
    private static final String DEFAULT_BATCH_COL = "NULL";

    private int randomizations = DEFAULT_RANDOMIZATIONS;
    private String posteriors = DEFAULT_POSTERIORS;
    private String batchCol = DEFAULT_BATCH_COL;

    //
    // Getters
    //

    /**
     * Get the number of randomizations to calculate
     *
     * @return Returns the number of randomizations (int)
     */
    protected int getRandomizations() {
        return this.randomizations;
    }

    /**
     * Get option for model plotting
     *
     * @return option for model plotting
     */
    protected String getModelPlotOption() {
        return this.modelPlot;
    }

    /**
     * Get Docker image to run
     *
     * @return Returns docker image (String)
     */
    protected String getDockerImage() {
        return this.docker;
    }

    /**
     * Get posteriors option
     *
     * @return Returns posteriors option (String)
     */
    protected String getPosteriorsOption() {
        return this.posteriors;
    }

    /**
     * Get group option for model fitting
     *
     * @return Returns group column name (String)
     */
    protected String getGroupCol() {
        return this.groupCol;
    }

    /**
     * Get minimum observations for using one gene during model fitting
     *
     * @return Returns minimum observation to keep a gene (int)
     */
    protected int getMinObservations() {
        return this.minObservation;
    }

    /**
     * Get minimum number of genes to use for model fitting
     *
     * @return Returns minimum number of genes for model fitting (int)
     */
    protected int getMinGenes() {
        return this.minGenes;
    }

    /**
     * Get option for quick threshold segmentation of dropped out features
     *
     * @return option for threshold segmentation (String)
     */
    protected String getThresholdSegOption() {
        return this.thresholdSeg;
    }

    /**
     * Get threshold for identification of dropped out features
     *
     * @return threshold for segmentation (int)
     */
    protected int getFailureThreshold() {
        return this.failureThreshold;
    }

    /**
     * Get max pairs for comparison during model fitting
     *
     * @return maximum number of pairs for model fitting (int)
     */
    protected int getMaxPairs() {
        return this.maxPairs;
    }

    /**
     * Get min pairs for comparison during model fitting
     *
     * @return minimum number of pairs for model fitting (int)
     */
    protected int getMinPairs() {
        return this.minPairs;
    }

    /**
     * Get linear fit option for dispersion estimate
     *
     * @return option for linear fitting of dispersion parameter (String)
     */
    protected String getLinearFitOption() {
        return this.linearFit;
    }

    /**
     * Get poisson distribution parameter (for drop-out component)
     *
     * @return poisson parameter (double)
     */
    protected double getPoissonPar() {
        return this.poissonParam;
    }

    /**
     * Get minimum value for dispersion parameter
     *
     * @return minimum value of dispersion parameter (double)
     */
    protected double getMinTheta() {
        return this.minTheta;
    }

    /**
     * Get maximum value for dispersion parameter
     *
     * @return maximum value of dispersion parameter (double)
     */
    protected double getMaxTheta() {
        return this.maxTheta;
    }

    /**
     * Get number of points for prior distribution estimate
     *
     * @return number of points for prior estimate (int)
     */
    protected int getPriorLength() {
        return this.length;
    }

    /**
     * Get prior plot option
     *
     * @return prior plot option (String)
     */
    protected String getPriorPlotOption() {
        return this.priorPlot;
    }

    /**
     * Get pseudocount for prior estimate
     *
     * @return pseudocount for prior estimate (int)
     */
    protected int getPseudocount() {
        return this.pseudocount;
    }

    /**
     * Get quantile to use for maximum expression estimate
     *
     * @return quantile for maximum expression estimate (double)
     */
    protected double getQuantile() {
        return this.quantile;
    }

    /**
     * Get value for maximum expression estimate
     *
     * @return maximum value for expression estimate (Integer)
     */
    protected Integer getMaxValue() {
        return this.maxValue;
    }

    /**
     * Get column name for batch conditions
     *
     * @return column name for batch conditions (String)
     */
    protected String getBatchCol() {
        return this.batchCol;
    }

    //
    // Module Methods
    //
    @Override public InputPorts getInputPorts() {
        return new InputPortsBuilder().addPort("sce", false, this.SCE_FORMAT).create();
    }

    @Override public String getName() {
        return MODULE_NAME;
    }

    @Override public Version getVersion() {
        return Globals.APP_VERSION;
    }

    @Override public OutputPorts getOutputPorts() {
        return singleOutputPort("dgeoutput", RESULT_FORMAT);
    }

    @Override public void configure(final StepConfigurationContext context,
                                    final Set<Parameter> stepParameters) throws EoulsanException {

        for (Parameter p : stepParameters) {

            switch (p.getName()) {

                case DOCKER_IMAGE:
                    this.docker = p.getStringValue();
                    break;

                case CORES:
                    this.nCores = p.getIntValue();
                    break;

                case MODEL_FIT_COL:
                    this.groupCol = p.getStringValue();
                    break;

                case PRIOR_LENGTH:
                    this.length = p.getIntValue();
                    break;

                case BATCH_COL:
                    this.batchCol = p.getStringValue();
                    break;

                case RANDOMIZATIONS:
                    this.randomizations = p.getIntValue();
                    break;

                case MIN_OBS:
                    this.minObservation = p.getIntValue();
                    break;

                case MIN_GENES:
                    this.minGenes = p.getIntValue();
                    break;

                case THRESHOLD_SEG:
                    this.thresholdSeg =
                            Boolean.toString(p.getBooleanValue()).toUpperCase();
                    break;

                case FAILURE_THRESHOLD:
                    this.failureThreshold = p.getIntValue();
                    break;

                case MAX_PAIRS:
                    this.maxPairs = p.getIntValue();
                    break;

                case MIN_PAIRS:
                    this.minPairs = p.getIntValue();
                    break;

                case POISSON_PARAM:
                    this.poissonParam = p.getDoubleValue();
                    break;

                case LINEAR_FIT:
                    this.linearFit =
                            Boolean.toString(p.getBooleanValue()).toUpperCase();
                    break;

                case MIN_THETA:
                    this.minTheta = p.getDoubleValue();
                    break;

                case MAX_THETA:
                    this.maxTheta = p.getDoubleValue();
                    break;

                case MODEL_PLOTS:
                    this.modelPlot = Boolean.toString(p.getBooleanValue());
                    break;

                //case CROSSFIT_PLOTS:
                //this.crossfitPlot = Boolean.toString(p.getBooleanValue());
                //break;

                case PRIOR_PLOT:
                    this.priorPlot = Boolean.toString(p.getBooleanValue());
                    break;

                case PSEUDOCOUNT:
                    this.pseudocount = p.getIntValue();
                    break;

                case QUANTILE:
                    this.quantile = p.getDoubleValue();
                    break;

                case MAX_VALUE:
                    this.maxValue = p.getIntValue();
                    break;

                case POSTERIORS:
                    this.priorPlot = Boolean.toString(p.getBooleanValue());
                    break;

                default:
                    Modules.unknownParameter(context, p);
            }
        }

        // Log Step parameters
        getLogger().info(
                "In " + getName() + ", general options: " + "image=" + this.docker
                        + ", n.cores=" + this.nCores);

        getLogger().info("In " + getName() + ", for model fitting: "
                + "\n data manipulation: " + "model.group.col=" + this.groupCol
                + ", min.observation=" + this.minObservation + ", min.genes="
                + this.minGenes + "\n dropout estimation: "
                + ", threshold.segmentation=" + this.thresholdSeg
                + ", failure.threshold=" + this.failureThreshold + ", max.pairs="
                + this.maxPairs + ", min.pairs=" + this.minPairs
                + ", poisson.param=" + this.poissonParam
                + "\n Negative binomial estimation: " + "linear.fit="
                + this.linearFit + ", min.theta=" + this.minTheta + ", max.theta="
                + this.maxTheta + "\n save.model.plot=" + this.modelPlot);

        getLogger().info(
                "In " + getName() + ", for prior calculation: " + "prior.length="
                        + this.length + ", pseudocount=" + this.pseudocount
                        + ", quantile=" + this.quantile + ", max.value=" + this.maxValue
                        + "save.prior.plot= " + this.priorPlot);

        getLogger().info(
                "In " + getName() + ", for test: " + "batch.col=" + this.batchCol
                        + ", n.randomizations=" + this.randomizations
                        + ", return.posteriors" + this.posteriors);

        // Set executor
        this.executor = new DockerExecutorInterpreter(docker);
    }

    @Override public TaskResult execute(final TaskContext context,
                                        final TaskStatus status) {

        // Get Design
        final Design design = context.getWorkflow().getDesign();

        // Get Input files
        final String name = context.getInputData(SCE_FORMAT).getName();
        final String format =
                context.getInputData(SCE_FORMAT).getDataFile().toFile().getPath();
        getLogger().info(
                "Input: " + format + " " + name);
        final String sce =
                context.getInputData(SCE_FORMAT).getDataFile().toFile().getPath();
        //final String InputGenes = context.getInputData(GENES_FORMAT).getDataFile().toFile().getPath();

        // Construct Prefix
        final String StepPrefix = context.getCurrentStep().getId() + "_";
        final String OutputModel = StepPrefix + MODEL_PREFIX + ".tsv";
        final String OutputPrior = StepPrefix + PRIOR_PREFIX + ".rds";

        // Get Experiments
        List<Experiment> experiments = design.getExperiments();

        //Set Output Files
        File execDir = context.getStepOutputDirectory().toFile();
        File OutputDir = context.getOutputDirectory().toFile();
        File logDirectory =
                ((TaskContextImpl) context).getTaskOutputDirectory().toFile();
        File tmpDir = context.getLocalTempDirectory();
        File stdoutFile = new File(logDirectory, "scde" + ".model.out");
        File stderrFile = new File(logDirectory, "scde" + ".model.err");

        try {

            // Calculate error models for cells see toolExecutor for galaxytools
            String commandLine =
                    buildFittingCommandLine(sce, OutputModel);
            List<String> command = executor.createCommandLine(commandLine);
            executor.execute(command, execDir, tmpDir, stdoutFile, stderrFile,
                    new File[] {OutputDir});

            // Calculate prior distribution
            stdoutFile = new File(logDirectory, "scde" + ".prior.out");
            stderrFile = new File(logDirectory, "scde" + ".prior.err");

            commandLine = buildPriorCommandLine(OutputModel, sce,
                    OutputPrior);
            command = executor.createCommandLine(commandLine);
            executor.execute(command, execDir, tmpDir, stdoutFile, stderrFile,
                    new File[] {OutputDir});

            // Launch corresponding function : list of compa or classic (all to ref)
            for (Experiment e : experiments) {

                // Get Experiment name, columns (model) and comparisons
                String expID = e.getId();
                String expName = e.getName();
                String columns;
                ExperimentMetadata md = e.getMetadata();

                if (md.contains("columns")) {
                    columns = md.get("columns");
                } else {
                    columns = e.getMetadata().getModel();
                    columns = columns.replace("~", "");
                    getLogger().warning("In " + getName() + " Exp." + expName
                            + ": columns attribute not found, used formula instead.");
                }

                List<String> comparisons =
                        Arrays.asList(e.getMetadata().getComparisons().split(";"));

                // For each comparison
                for (String comparison : comparisons) {

                    // Split name and content
                    String[] table = comparison.split(":");

                    // Check if comparison is well formed
                    if (table.length > 2)
                        throw new EoulsanException("Invalid comparison format");

                    // Get separated name and content
                    String compName = table[0];
                    String c = table[1];

                    // Create output name
                    String OutputResult =
                            StepPrefix + expName + "_" + compName + ".tsv";

                    // Run test
                    stdoutFile =
                            new File(logDirectory, "scde" + "." + expID + ".out");
                    stderrFile =
                            new File(logDirectory, "scde" + "." + expID + ".err");
                    commandLine = buildTestCommandLine(OutputModel, sce, OutputPrior, columns, expID, c,
                            OutputResult);
                    command = executor.createCommandLine(commandLine);
                    executor.execute(command, execDir, tmpDir, stdoutFile,
                            stderrFile, new File[] {OutputDir});
                }
            }

        } catch (Exception e) {
            return status.createTaskResult(e,
                    "Error while running SCDE : " + e.getMessage());
        }
        return status.createTaskResult();
    }

    /**
     * Command builder for model fitting script
     *
     * @param sce: path to count Single Cell Experiment (String)
     * @param modelFile:  path to output model file (String)
     * @return command line (String)
     */
    protected final String buildFittingCommandLine(String sce, String modelFile) {
        return (cmd + " " + modelScript + " " + sce + " "
                + Integer.toString(getNbCores()) + " " + getGroupCol() + " "
                + getModelPlotOption() + " " + Integer
                .toString(getMinObservations()) + " " + Integer
                .toString(getMinGenes()) + " " + getThresholdSegOption() + " "
                + Integer.toString(getFailureThreshold()) + " " + Integer
                .toString(getMaxPairs()) + " " + Integer.toString(getMinPairs())
                + " " + Double.toString(getPoissonPar()) + " "
                + getLinearFitOption() + " " + Double.toString(getMinTheta()) + " "
                + Double.toString(getMaxTheta()) + " " + modelFile);
    }

    /**
     * Command builder for prior script
     *
     * @param modelFile:  path to model file (String)
     * @param sce: path to count Single Cell Experiment (String)
     * @param priorFile:  path to output prior file (String)
     * @return command line (String)
     */
    protected final String buildPriorCommandLine(String modelFile,
                                                 String sce, String priorFile) {
        return (cmd + " " + priorScript + " " + modelFile + " " + sce
                + " " + Integer.toString(getPriorLength()) + " "
                + getPriorPlotOption() + " " + Integer.toString(getPseudocount())
                + " " + (getQuantile() == 0 ? "NULL" :
                Double.toString(getQuantile())) + " " + (getMaxValue() != null ?
                getMaxValue().toString() : "NULL") + " " + priorFile);
    }

    /**
     * Command builder for test script
     *
     * @param modelFile:  path to model file (String)
     * @param sce: path to count Single Cell Experiment (String)
     * @param priorFile:  path to prior file (String)
     * @param columns:    columns for comparison (String)
     * @param expID:      ID of experience (String)
     * @param comparison: comparison (String)
     * @param outputFile: path to result file (String)
     * @return command line (String)
     */
    protected final String buildTestCommandLine(String modelFile,
                                                String sce, String priorFile, String columns,
                                                String expID, String comparison, String outputFile) {
        return (cmd + " " + testScript + " " + modelFile + " " + sce + " " + priorFile + " " + columns + " " + expID
                + " " + comparison + " " + Integer.toString(this.getNbCores()) + " "
                + Integer.toString(this.getRandomizations()) + " " + getBatchCol()
                + " " + this.getPosteriorsOption() + " " + outputFile);
    }
}
