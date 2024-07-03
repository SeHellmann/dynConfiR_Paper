Contributor(s): 
- LAW, Frankie Ho Fai
- LEE, Alan Lap Fai (ALANLEE@LN.EDU.HK)

Citation: unpublished, but presented in a talk at ASSC 23

<Details of experiment>

Stimulus: random-dot kinematograms generated using PsychoPy2

Task: simultaneous direction-discrimination (left or right) and 4-point, confidence-rating

Confidence scale: 
- a 4-point confidence rating, simultaneously indicated with perceptual decisions (left-right direction discrimination) on a US-English keyboard
- [z, x, c, v] to indicate LEFTWARD motion [z: "Surely LEFT", v: "Guess LEFT"]
- [m, comma, period, slash] to indicate RIGHTWARD motion [m: "Guess RIGHT", slash: "Surely RIGHT"]

Manipulation: 
- task difficulty controlled by varying motion coherence (noise dots moved in random directions)

Calibration (first 240 trials):
-  randomly interleaved with 6 independent staircases (40 trials each; algorithm: ASA by Kesten, 1958)
- 3 targeted levels of direction-discrimination accuracy (1: .52, 3: .65, and 5: .78; see accTarget in data)
- each targeted level had 2 independent staircases (see stairID in data)
- coherence thresholds for levels 1, 3, and 5 were then estimated by averaging the last trial's estimates of the two staircases;
- threshold for coherence level 2 was estimated by averaging those between 1 and 3
- threshold for coherence level 4 was estimated by averaging those between 3 and 5

Block size:
- 60 trials per block
- 20 blocks per subject
- subjects took a brief break after each block

Feedback: 
- no trial-by-trial feedback
- block feedback was given after each block, which includes two accuracies (both as percentages of correct trials) :
	1. accuracy in the immediately-past block
	2. the cumulative accuracy since the beginning of the experiment

<Key to variables in data>
Subj_idx        : subject ID (1 - 16)
Stimulus        : direction of RDK, 1:left, 2:right
Response        : perceptual response, 1:left, 2:right
Confidence      : confidence rating, 1:lowest, 4:highest
RT_decConf      : response time (in seconds) of the simultaneous decision-and-confidence response
blocki          : the ith block in the experiment (60 trials per block)
triali          : the ith trial in the experiment (1200 trials per subject)
keypress        : key pressed by the subject for this trial; see above for details
coh_level       : level of coherence, 1:lowest, 5: highest
coherence       : proportion of dots moving in the signal direction, 0: all dots moved randomly, 1:all dots moved coherently
isCalibTrial    : 1:calibration trial (i.e., first 240 trials), 0:not calibration; see above for details
stairID         : for calibration trial only (NaN for non-calibration trials); indicates which staircase for this coherence level
accTarget       : for calibration trial only (NaN for non-calibration trials); targeted accuracy level for this staircase