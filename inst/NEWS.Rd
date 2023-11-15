\name{NEWS}
\title{News for \R Package \pkg{Pigenegene}}

\section{Changes in version 1.29.8 (2023-11-13)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item Using "/" is avoided in  \code{determine.modules} and
      \code{learn.bn} functions to support running them on Windows.
      The \code{"toLocal"} task in \code{learn.bn()} and \code{sbatch()}
      still do not support Windows.
    }
  }
  \subsection{General}{
    \itemize{
      \item \code{inst/script/bn.calculation.job.R} is added as an
      example for the \code{bnCalculationJob} argument. 
    }
  }
}

\section{Changes in version 1.29.6 (2023-11-10)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item Setting \code{toCompact} to \code{TRUE} in
      \code{make.decision.tree()} has the same effect as \code{NULL}.
    }
  }
}

\section{Changes in version 1.29.4 (2023-11-09)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item The \code{averaged.network()} function now correctly
      passed test data to the \code{module.heatmap()} function. 
    }
  }
  \subsection{Changes in existing functions}{
    \itemize{
      \item More details on the effect of test data in the documentation of
      \code{one.step.pigengene()} and \code{make.decision.tree()}
      functions. \code{make.decision.tree()} returns more information on
      performance.
    }
  }
}

\section{Changes in version 1.27.16 (2023-06-16)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item Scaling is done in the \code{compute.pigengene()} function
      to avoid the following error that Sogand detailes on her lano:
      Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'
    }
  }
}

\section{Changes in version 1.27.2 (2023-05-12)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The \code{compute.pigengene()} function does not warn about
      almost constant genes.
    }
  }
}

\section{Changes in version 1.25.16 (2023-03-22)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item Constant genes are now ignored in the
      \code{compute.pigengene(doWgcna=FALSE,...)} function to prevent a
      run time error. 
    }
  }
}

\section{Changes in version 1.25.12 (2023-03-02)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The default value changed in
      \code{gene.mapping(leaveNA=FALSE, ...)}.
    }
  }
}

\section{Changes in version 1.25.10 (2023-01-04)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The \code{doWgcna} option is added to the
      \code{compute.pigengene} function.
    }
  }
}

\section{Changes in version 1.25.4 (2022-12-01)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The \code{get.enriched.pw} function
      adds gene symbols in the excel file.
    }
  }
}

\section{Changes in version 1.23.4 (2022-05-23)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The beta table can now be suppressed in the
      \code{determine.modules} function.
    }
  }
}

\section{Changes in version 1.21.40 (2022-04-15)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item Habil changed the default value from
      \code{hu.mouse(host="useast.ensembl.org", ...)} to
      \code{hu.mouse(host="www.ensembl.org", ...)} to
      prevent a check error on Bioconductor.
    }
  }
}
  
\section{Changes in version 1.21.36 (2021-11-16)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item Habil added the \code{doReturNetworks} argument to \code{one.step.pigengene()}.
    }
  }
}

\section{Changes in version 1.21.34 (2021-11-12)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item Habil renamed \code{identify.modules()} to
      \code{determine.modules()}.
    }
  }
}

\section{Changes in version 1.21.30 (2021-11-12)}{
  \subsection{New functions}{
    \itemize{
      \item Neda exported \code{identify.modules()},
      \code{make.filter()}, and \code{apply.filter()} functions.
    }
  }
}

\section{Changes in version 1.19.50 (2021-10-14)}{
  \subsection{New functions}{
    \itemize{
      \item Neda added the \code{identify.modules()},
      \code{make.filter()}, and \code{apply.filter()} functions, but not
      exported them yet. 
    }
  }
}

\section{Changes in version 1.19.30 (2021-09-08)}{
   \subsection{Bug Fixes}{
       \itemize{
         \item The \code{averaged.network()} function does not
	 have the \code{nodes} argument in bnlearn Version >=4.7, and
	 thus this argument was removed.
    }
  }
}

\section{Changes in version 1.19.24 (2021-08-06)}{
  \subsection{New functions}{
    \itemize{
      \item The \code{get.enriched.pw()} function is added. 
    }
  }
   \subsection{Bug Fixes}{
       \itemize{
         \item A bug fix in the \code{gene.mapping()} function that used to occur
         when we had multiple output databases. 
    }
  }
}

\section{Changes in version 1.19.10 (2021-06-25)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The \code{message.if()} function can now write the message
      in a text file.
    }
  }
}

\section{Changes in version 1.19.8 (2021-05-25)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item The C50 plot function seems to have different behavior when
      the number of Labels is 2. Habil reverse the color to fix the
      resulting bug.
    }
  }
}

\section{Changes in version 1.17.10 (2021-02-02)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item When building the Ubuntu, Neda got an error message on
      building the vignette "argument is of length zero". It was due to
      her old version of BiocStyle, 2.14.4. To solve this issue,
      Pigengene now depends BiocStyle Version >= 2.19.1. Even 2.18.1
      might be enough with R 4.0.3.
    }
  }
}

\section{Changes in version 1.15.20 (2020-09-28)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item When automatically computing the \code{threshold} in the
      \code{consensus} function, we now do not allow it to be more than
      1 to be compatible with bnlearn 4.6.1.
    }
  }
}

\section{Changes in version 1.15.16 (2020-09-22)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The modules can now be determined in the
      \code{module.heatmap} function using the new \code{mes} argument.
    }
  }
}

\section{Changes in version 1.15.14 (2020-09-08)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{combine.networks} now really removes the big TOM file
      when \code{doRemoveTOM=TRUE}.
    }
  }
}

\section{Changes in version 1.15.12 (2020-06-22)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The \code{DiseaseColors} argument in the
      \code{plot.pigengene} function can now be set automatically.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{pheatmap.type} now works fine even when the number of
      samples in a condition is only 1.
    }
  }
}

\section{Changes in version 1.13.44 (2020-04-21)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item \code{one.step.pigengene} now talks more, and automatically
      converts input data frames to matrices.
      \item \code{module.heatmap} now talks more, and does not get stuck
      if a module has only one gene, or if the expression levels of some
      genes in a module are not available in the data.
    }
  }
}

\section{Changes in version 1.13.14 (2020-03-03)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item \code{inherits(a, b)} is now used instead of \code{class(a)==b}.
    }
  }
}

\section{Changes in version 1.13.6 (2020-01-03)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item \code{doSave} was added to the \code{combine.networks()} function.
    }
  }
}

\section{Changes in version 1.13.4 (2019-11-19)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item In the \code{combine.networks()} function, the \code{netwok}
      typo was fixed.
    }
  }
}

\section{Changes in version 1.11.34 (2019-10-21)}{
  \subsection{General}{
    \itemize{
      \item The pipeline is now explained step by step in the vignette.  
    }
  }
}

\section{Changes in version 1.11.32 (2019-10-18)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The \code{Data} argument of \code{compute.pigengene} can now
      be a matrix with only 1 column. 
    }
  }
}

\section{Changes in version 1.11.30 (2019-10-02)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The \code{doRetuNetworks} argument is now added to the
      \code{combine.networks} function.
    }
  }
}

\section{Changes in version 1.11.28 (2019-10-01)}{
  \subsection{New functions}{
    \itemize{
      \item \code{message.if()} is now exported.
    }
  }
}

\section{Changes in version 1.11.26 (2019-09-26)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item Better QC in the \code{gene.mapping()} function, the
      possible keys will be printed if the input is not appropriate.
    }
  }
  \subsection{New functions}{
    \itemize{
      \item \code{save.if()} is now exported.
    }
  }
}

\section{Changes in version 1.11.24 (2019-09-03)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item In the \code{combine.networks()} function,
      \code{selectedModules} does not need to be \code{"All"}. Also, if
      \code{saveFile=NULL}, nothing will be saved without any error.
    }
  }
}

\section{Changes in version 1.11.20 (2019-05-15)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item \code{Data} and \code{Labels} can now be lists, which will be combined using 
      \code{combine.networks()} before analysis.
    }
  }
}

\section{Changes in version 1.11.4 (2019-05-02)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{repeat.data(times=1,...)} now produces valid output.
    }
  }
}

\section{Changes in version 1.9.26 (2019-04-24)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item In the \code{compute.pigengene} function, if a module has
      only one gene, its name is now not omitted in the csv file.
    }
  }
}

\section{Changes in version 1.9.24 (2019-04-24)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item An object of \code{\link{pigengene-class}} can now have a
      \code{heavyToLow} attribute.
    }
  }
}

\section{Changes in version 1.9.20 (2019-04-12)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item In the \code{compute.pigengene} function, the columns of
      \code{Data} and names of \code{modules} can now differ.
    }
  }
}

\section{Changes in version 1.9.14 (2019-02-12)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The \code{dOrderByW} argument is now added for the \code{compute.pigengene} function.
    }
  }
}

\section{Changes in version 1.9.8 (2018-11-22)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item All of the \code{learn.bn} function can now be set from
      \code{one.step.pigengene} through \code{bnArgs}.
    }
  }
}

\section{Changes in version 1.9.4 (2018-11-16)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item \code{gene.mapping} can now output multiple conventions.
    }
  }
}

\section{Changes in version 1.7.2 (2018-05-22)}{
  \subsection{General}{
    \itemize{
      \item The version of the package \code{C50} is now required to be
      at least 0.1.2, which exports the \code{as.party.C5.0()} function.
    }
  }
}

\section{Changes in version 1.5.22 (2018-04-28)}{
  \subsection{Bug Fixes}{
   \itemize{
      \item \code{stats::cor} is used in \code{compute.pigengene} and
  \code{draw.cor} functions. See the NAMESPACE for the important reason.
    }
  }
}

\section{Changes in version 1.5.9 (2018-03-12)}{
  \subsection{General}{
    \itemize{
      \item \code{cor} is imported from WGCNA, but not from stats, because WGCNA does
      not call the \code{cor} function properly.
    }
  }
  \subsection{Changes in existing functions}{
    \itemize{
      \item \code{RsquaredCut} is added to the arguments of the \code{one.step.pigengene} function.
    }
  }
}

\section{Changes in version 1.5.6 (2018-01-19)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item In the \code{get.fitted.leaf} function, the function
  C50:::as.party.C5.0 is used, which used to be exported in the previous
  versions of C50, but not in version 0.1.1.
    }
  }
}

\section{Changes in version 1.5.2 (2017-11-10)}{
  \subsection{New functions}{
    \itemize{
      \item The \code{check.nas} function is now exported.
    }
  }
}

\section{Changes in version 1.3.8 (2017-09-13)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item Order of conditions in pheatmap.type can now be determined
  by the user.
    }
  }
}

\section{Changes in version 1.3.6 (2017-08-20)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item A bug in gene.mapping () function fixed to better map probe IDs.
    }
  }
}

\section{Changes in version 1.3.4 (2017-07-31)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The doTranspose argument added to the  heatmap.type() function.
    }
  }
}

\section{Changes in version 1.1.14 (2017-05-10)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item A bug in the module.heatmap() function fixed.
    }
  }
}

\section{Changes in version 1.1.12 (2017-04-15)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The mouseHomologFilter argument is added to the unexported
      furntion hu.mouse() with a default value of
      "with_mmusculus_homolog". This is now compatible with
      biomaRt_2.30.0, it used to be "with_homolog_mmus" in older versions of biomaRt.
    }
  }
}

\section{Changes in version 1.1.10 (2017-03-30)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item Checking the pigengene input of module.heatmap().
      \item Issues in the balance() function (not exported) where Labels
      is a factor are resolved. Also,  if all sampls have the same size,
      oversampling is automatically turned off.
      \item If Labels is a factor, it is now converted to a character
      vector in check. pigengene.input().
    }
  }
}

\section{Changes in version 1.1.6 (2017-03-27)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The module.heatmap() function now has the doAddEigengene
      and scalePngs arguments.
      \item The compute.pigengene()  function now reports also the size of
      modules in the pigengene_pvalue.csv output file.
    }
  }
}

\section{Changes in version 0.99.25 (2016-10-02)}{
  \subsection{Changes in existing functions}{
    \itemize{
      \item The compute.pigengene() function now uses welch.pvalue() instead of pvalues.manov().
    }
  }
}

\section{Changes in version 0.99.8 (2016-05-18, git f50187c)}{
  \subsection{General}{
    \itemize{
      \item  Under review by Bioconductor.
      \item Created.
    }
  }
}

%% \subsection{General}
%% \subsection{New functions}{
%% \subsection{Changes in existing functions} %% functionality, signature, output
%% \subsection{Bug Fixes}
