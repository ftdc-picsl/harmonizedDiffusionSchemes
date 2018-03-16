# Read INRIA samples
samples = read.csv("samples.csv")

# Check we have the expected number of directions on each shell
summary(factor(samples$shell))

# Source a bunch of R functions (can replace with package later)
functionFiles = list.files('/data/jet/pcook/code/diffusionGradientUtils', pattern = "*.R$", full.names = T)

for (f in functionFiles) {
  source(f)
}


inriaScheme = inriaShellsToBvecsBvals(samples, shells = c(1000,2000))

## Quick check
check = cbind(inriaScheme$bvecs, inriaScheme$bvals, samples)

## Check shells have the correct b-values
by(check[,4], factor(check[,5]), summary)


## Note that for the rest of the code, eg plotGradDirs, shells are assumed to include b=0, ie shells should be c(0,1000,2000) and b=0 is not used for plotting. We'll add some b=0 measurements later

## Look at b=1000 shell
plotGradDirs(inriaScheme$bvals, inriaScheme$bvecs, shells = c(0,1000,2000), whichShells = 2, opacity = 1, projectToUpper = T )

## Add b=2000 shell
plotGradDirs(inriaScheme$bvals, inriaScheme$bvecs, shells = c(0,1000,2000), whichShells = 3, opacity = 1, projectToUpper = T )

## Close the window and plot again without projection to the hemisphere
plotGradDirs(inriaScheme$bvals, inriaScheme$bvecs, shells = c(0,1000,2000), opacity = 1, projectToUpper = F )

## The INRIA dirs are spread on the full sphere, which is good

## Now read the 12-dir shell
elec12 = normalizePoints(read.csv("elec012.csv"))

## B=500 matches ADNI 3, though they only get 6
## Make bvecs a matrix for compatibility with other functions
elec12Scheme = list(bvecs = as.matrix(elec12), bvals = rep(500,12))

plotGradDirs(elec12Scheme$bvals, elec12Scheme$bvecs, shells = c(0,500), opacity = 1, projectToUpper = F)

## These directions are on a hemisphere, so let's distribute them
elec12VecsFullSphere = distributeShellsOverSphere(elec12Scheme$bvals, elec12Scheme$bvecs, shells = c(0,500))

## Check
cbind(elec12Scheme$bvecs, elec12VecsFullSphere)

## replace bvecs
elec12Scheme$bvecs = elec12VecsFullSphere

plotGradDirs(elec12Scheme$bvals, elec12Scheme$bvecs, shells = c(0,500), opacity = 1, projectToUpper = F)

## Now add some zeros
##
## HCP has N = 288 and adds M = 18 zeros, so M = N / 16, or 15.15 if you include the built-in 0
## at the start of the sequence
##
## M = 9 for 156 measurements is about the same
M = 9

zeroScheme = list(bvecs = matrix(rep(0,M*3), ncol = 3), bvals = rep(0,M))


## Merge the schemes
threeShellScheme = mergeSchemes(inriaScheme$bvals, inriaScheme$bvecs, elec12Scheme$bvals, elec12Scheme$bvecs)

threeShellsWithZeros = mergeSchemes(threeShellScheme$bvals, threeShellScheme$bvecs, zeroScheme$bvals, zeroScheme$bvecs, offsetFromEnd = F)

plotGradDirs(threeShellsWithZeros$bvals, threeShellsWithZeros$bvecs, shells = c(0,500,1000,2000), opacity = 1, projectToUpper = F)

scaled = writeSiemensGradientTable(threeShellsWithZeros$bvals, threeShellsWithZeros$bvecs, "siemens_M9_N156_multishell.txt")

by(scaled$gradientMagnitude, factor(threeShellsWithZeros$bvals), summary )



## Inner shell from the INRIA direction set
## Designed for fast DTI at b=1000

innerShellScheme = list(bvecs = inriaScheme$bvecs[which(inriaScheme$bvals == 1000),], bvals = inriaScheme$bvals[which(inriaScheme$bvals == 1000)])

## M = 4 for a total of 5 b=0 or a ratio of 9.6:1
M = 4

zeroScheme = list(bvecs = matrix(rep(0,M*3), ncol = 3), bvals = rep(0,M))

innerShellWithZeros = mergeSchemes(innerShellScheme$bvals, innerShellScheme$bvecs, zeroScheme$bvals, zeroScheme$bvecs, offsetFromEnd = F)

scaledInner = writeSiemensGradientTable(innerShellWithZeros$bvals, innerShellWithZeros$bvecs, "siemens_M4_N48_unishell.txt")

by(scaledInner$gradientMagnitude, factor(innerShellWithZeros$bvals), summary )


## Now produce an optimal 48-dir scheme for comparison

elec48 = normalizePoints(read.csv("elec048.csv"))

## B=500 matches ADNI 3, though they only get 6
## Make bvecs a matrix for compatibility with other functions
elec48Scheme = list(bvecs = as.matrix(elec48), bvals = rep(1000,48))

## These directions are on a hemisphere, so let's distribute them
elec48VecsFullSphere = distributeShellsOverSphere(elec48Scheme$bvals, elec48Scheme$bvecs, shells = c(0,1000))

## Check
cbind(elec48Scheme$bvecs, elec48VecsFullSphere)

## replace bvecs
elec48Scheme$bvecs = elec48VecsFullSphere

plotGradDirs(elec48Scheme$bvals, elec48Scheme$bvecs, shells = c(0,1000), opacity = 1, projectToUpper = F)

elec48WithZeros = mergeSchemes(elec48Scheme$bvals, elec48Scheme$bvecs, zeroScheme$bvals, zeroScheme$bvecs, offsetFromEnd = F)

scaledElec48 = writeSiemensGradientTable(elec48WithZeros$bvals, elec48WithZeros$bvecs, "siemens_M4_N48_elec.txt")

by(scaledElec48$gradientMagnitude, factor(elec48WithZeros$bvals), summary )
