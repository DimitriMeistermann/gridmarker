#' Project features linked to a projection from each point on a grid
#'
#' @param proj Matrix or dataframe containing at least two vectors of numeric corresponding to the x/y coordinates.
#' @param featureMat A matrix with row as features (if `transpose`).
#' @param gridResolution Increase or decrease the number of points in the grid. It corresponds to the number of grid points along one axis if the projection is a square.
#' @param method Method for computing the feature specificity score. Can be `"aurocFC"`  (auroc * fold-change), `"auroc"` or `"FC"` (Fold-change).
#' @param grid A dataframe returned by this function if `returnGrid`. Contains data for plotting the grid (see `returnGrid` for more explanations)
#' @param transpose If `transpose`, samples are columns and features are rows.
#' @param geomTextFun Function used to plot the text if returnGrid=FALSE.
#' @param fontface Font face of the displayed features. Passed to the function specified in `geomTextFun`
#' @param returnPseudoSamples Return the convoluted samples that correspond to each point of the grid.
#' @param returnGrid
#' Return a dataframe instead of plotting, corresponding to the coordinate of the grid,
#' Contains the following columns: `X` and `Y` (coordinates), `nContrib` (number of contributing observation to the grid point, the grid point is not created if `nContrib==0`),
#' `best` (best feature of the grid point), `score` (specificity score of the best feature), `sizeScore` (score weighted by number of observations, will determine the fontsize).
#' @param ... additionnal parameters passed to the text plotting function-
#'
#' @details
#' 1. Project the dataset on a 2D plane using a suitable algorithm (e.g. PCA or t-SNE).
#' 2. Create a grid of points with the specified grid_size.
#' 3. For each point on the grid, find its neighborhood by looking at the points in the dataset that are closest to it.
#' 4. For each feature of the featureMat, calculate the specificity of that feature in the neighborhood of the grid point. This could be done using a metric such as the average distance between points that have that feature and points that do not have that feature.
#' 5. Display the feature name with the highest specificity for each grid point.
#'
#' @return A ggproto object to add to a 2d projection by default. Can also be a dataframe if `returnGrid` or `returnPseudoSamples`.
#' @export
#'
#' @examples
#' library(ggplot2)
#' library(oob)
#' data(iris)
#'
#' qDataIris<-t(iris[,1:4])
#' irisUMAP<-UMAP(qDataIris)
#' ggproj<-ggplot(data.frame(x=irisUMAP[,1],y=irisUMAP[,2]),aes(x=x,y=y))+geom_point()+
#'   theme(panel.background = element_rect(fill = NA,
#'     colour = "black"), panel.grid.major = element_line(colour = NA),
#'     panel.grid.minor = element_line(colour = NA)) + guides(size = "none")
#' ggproj+markergrid(irisUMAP, qDataIris)+scale_size_continuous(range = c(1,3))
#' ggproj+markergrid(irisUMAP, qDataIris, geomTextFun=geom_text)+scale_size_continuous(range = c(1,3))
#' grid<-markergrid(irisUMAP,qDataIris, returnGrid=TRUE)
markergrid <-
    function(proj,
            featureMat,
            gridResolution = 15,
            method = c("aurocFC", "auroc", "FC"),
            grid = NULL,
            transpose = TRUE,
            geomTextFun = shadowtext::geom_shadowtext,
            fontface = "bold.italic",
            returnPseudoSamples = FALSE,
            returnGrid = FALSE,
            ...) {
        if (transpose)
            featureMat <- t(featureMat)
        if (is.null(colnames(featureMat)))
            stop("features should be named")
        method <- method[1]
        if (!method %in% c("aurocFC", "auroc", "FC"))
            stop("method should have one of the following value: 'aurocFC', 'auroc' or 'FC'")
        if (is.null(grid)) {
            grid<-buildEmptyGrid(proj, gridResolution)

            pseudoSamples <-
                matrix(
                    data = 0,
                    nrow = nrow(grid),
                    ncol = ncol(featureMat)
                )
            colnames(pseudoSamples) <- colnames(featureMat)

            hasVal <- foreach(i = seq_along(nrow(grid))) %do% {
                centerX <- grid[i, 1]
                centerY <- grid[i, 2]

                isContrib2Grid <-
                    (proj[, 1] - centerX) ^ 2 + (proj[, 2] - centerY) ^ 2 < radius ^ 2
                grid[i, "nContrib"] <- sum(isContrib2Grid)

                if (sum(isContrib2Grid) == 0)
                    return(FALSE)

                isContribX <- proj[isContrib2Grid, 1]
                isContribY <- proj[isContrib2Grid, 2]

                dist2center <-
                    sqrt((isContribX - centerX) ^ 2 + (isContribY - centerY) ^ 2)

                w <-  1 - (dist2center / radius)
                pseudoSamples[i, ] <-
                    apply(featureMat[isContrib2Grid, , drop = FALSE], 2, function(x)
                        stats::weighted.mean(x, w))

                return(TRUE)
            }
            hasVal <- unlist(hasVal)

            grid <- grid[hasVal, ]
            pseudoSamples <- pseudoSamples[hasVal, , drop = FALSE]

            if (returnPseudoSamples)
                return(pseudoSamples)

            pseudoSamples <- scale(pseudoSamples, center = TRUE, scale = FALSE)
            rowIndex <- seq_along(nrow(pseudoSamples))

            if (method == "aurocFC" | method == "auroc")
                aurocs <- apply(pseudoSamples, 2, function(f) {
                    scoreRank <- rank(f)
                    n1 <- length(f) - 1

                    res <-
                        vapply(seq_along(length(scoreRank)), FUN.VALUE = numeric(1), function(i)
                            sum(scoreRank[-i]))
                    1 - (res - n1 * (n1 + 1) / 2) / n1
                })

            if (method == "aurocFC") {
                pseudoSamples <- pseudoSamples * aurocs
            } else if (method == "auroc") {
                pseudoSamples <- aurocs
            }

            maxValIndex <- apply(pseudoSamples, 1, which.max)

            grid[, "best"] <- colnames(pseudoSamples)[maxValIndex]

            grid$score <-
                NULL
            for (i in seq_along(nrow(pseudoSamples)))
                grid[i, "score"] <-  pseudoSamples[i, maxValIndex[i]]
            grid$sizeScore <-  grid$score * grid$nContrib
        }


        if (returnGrid)
            return(grid)

        return(geomTextFun(
            data = grid,
            inherit.aes = FALSE,
            mapping = aes_string(
                x = "X",
                y = "Y",
                label = "best",
                size = "sizeScore"
            ),
            fontface = fontface,
            ...
        ))
    }

buildEmptyGrid<-function(proj, gridResolution){
    minX <- min(proj[, 1])
    maxX <- max(proj[, 1])
    minY <- min(proj[, 2])
    maxY <- max(proj[, 2])

    xrange <- maxX - minX
    yrange <- maxY - minY

    radius <- sqrt(xrange * yrange) / gridResolution

    gridX <- seq(minX - radius, maxX + radius, radius)
    gridY <- seq(minY - radius, maxY + radius, radius)

    grid <- data.frame(expand.grid(gridX, gridY))
    colnames(grid) <- c("X", "Y")

    grid$nContrib <- NULL
    return(grid)
}
