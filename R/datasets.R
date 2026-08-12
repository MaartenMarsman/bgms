#' PTSD Symptoms in Wenchuan Earthquake Survivors Who Lost a Child
#'
#' This dataset contains responses to 17 items assessing symptoms of post-traumatic stress disorder (PTSD)
#' in Chinese adults who survived the 2008 Wenchuan earthquake and lost at least one child in the disaster
#' \insertCite{McNallyEtAl_2015}{bgms}. Participants completed the civilian version of the Posttraumatic Checklist,
#' with each item corresponding to a DSM-IV PTSD symptom. Items were rated on a 5-point Likert scale from
#' "not at all" to "extremely," indicating the degree to which the symptom bothered the respondent in the
#' past month.
#'
#' @format A matrix with 362 rows and 17 columns. Each row represents a participant.
#' \describe{
#'   \item{intrusion}{Repeated, disturbing memories, thoughts, or images of a
#'   stressful experience from the past?}
#'   \item{dreams}{Repeated, disturbing dreams of a stressful experience from
#'   the past?}
#'   \item{flash}{Suddenly acting or feeling as if a stressful experience were
#'   happening again (as if you were reliving it)?}
#'   \item{upset}{Feeling very upset when something reminded you of a stressful
#'   experience from the past?}
#'   \item{physior}{Having physical reactions (e.g., heart pounding, trouble
#'   breathing, sweating) when something reminded you of a stressful experience
#'   from the past?}
#'   \item{avoidth}{Avoiding thinking about or talking about a stressful
#'   experience from the past or avoiding having feelings related to it?}
#'   \item{avoidact}{Avoiding activities or situations because they reminded you
#'   of a stressful experience from the past?}
#'   \item{amnesia}{Trouble remembering important parts of a stressful
#'   experience from the past?}
#'   \item{lossint}{Loss of interest in activities that you used to enjoy?}
#'   \item{distant}{Feeling distant or cut off from other people?}
#'   \item{numb}{Feeling emotionally numb or being unable to have loving
#'   feelings for those close to you?}
#'   \item{future}{Feeling as if your future will somehow be cut short?}
#'   \item{sleep}{Trouble falling or staying asleep?}
#'   \item{anger}{Feeling irritable or having angry outbursts?}
#'   \item{concen}{Having difficulty concentrating?}
#'   \item{hyper}{Being "super-alert" or watchful or on guard?}
#'   \item{startle}{Feeling jumpy or easily startled?}
#' }
#'
#' @source \url{https://psychosystems.org/wp-content/uploads/2014/10/Wenchuan.csv}
#'
#' @docType data
#' @keywords datasets
#' @name Wenchuan
#' @usage data("Wenchuan")
#' @references
#' \insertAllCited{}
NULL

#' ADHD Symptom Checklist for Children Aged 6--8 Years
#'
#' This dataset includes ADHD symptom ratings for 355 children aged 6 to 8 years from the
#' Children's Attention Project (CAP) cohort \insertCite{Silk_2019_ADHD}{bgms}. The sample
#' consists of 146 children diagnosed with ADHD and 209 without a diagnosis. Symptoms were
#' assessed through structured interviews with parents using the NIMH Diagnostic Interview
#' Schedule for Children IV (DISC-IV) \insertCite{Shaffer_2000_nimh}{bgms}. The checklist
#' includes 18 items: 9 Inattentive (I) and 9 Hyperactive/Impulsive (HI). Each item is binary
#' (1 = present, 0 = absent).
#'
#' @format A data frame with 355 rows and 19 columns.
#' \describe{
#'   \item{group}{ADHD diagnosis: 1 = diagnosed, 0 = not diagnosed}
#'   \item{avoid}{Often avoids, dislikes, or is reluctant to engage in tasks
#'   that require sustained mental effort (I)}
#'   \item{closeatt}{Often fails to give close attention to details or makes
#'   careless mistakes in schoolwork, work, or other activities (I)}
#'   \item{distract}{Is often easily distracted by extraneous stimuli (I)}
#'   \item{forget}{Is often forgetful in daily activities (I)}
#'   \item{instruct}{Often does not follow through on instructions and fails to
#'   finish schoolwork, chores, or duties in the workplace (I)}
#'   \item{listen}{Often does not seem to listen when spoken to directly
#'   (I)}
#'   \item{loses}{Often loses things necessary for tasks or activities (I)}
#'   \item{org}{Often has difficulty organizing tasks and activities (I)}
#'   \item{susatt}{Often has difficulty sustaining attention in tasks or play
#'   activities (I)}
#'   \item{blurts}{Often blurts out answers before questions have been completed
#'   (HI)}
#'   \item{fidget}{Often fidgets with hands or feet or squirms in seat
#'   (HI)}
#'   \item{interrupt}{Often interrupts or intrudes on others (HI)}
#'   \item{motor}{Is often "on the go" or often acts as if "driven by a motor"
#'   (HI)}
#'   \item{quiet}{Often has difficulty playing or engaging in leisure activities
#'   quietly (HI)}
#'   \item{runs}{Often runs about or climbs excessively in situations in which
#'   it is inappropriate (HI)}
#'   \item{seat}{Often leaves seat in classroom or in other situations in which
#'   remaining seated is expected (HI)}
#'   \item{talks}{Often talks excessively (HI)}
#'   \item{turn}{Often has difficulty awaiting turn (HI)}
#' }
#'
#' @source \insertCite{Silk_2019_ADHD;textual}{bgms}.
#' Data retrieved from \doi{10.1371/journal.pone.0211053.s004}.
#' Licensed under the CC-BY 4.0: \url{https://creativecommons.org/licenses/by/4.0/}
#'
#' @docType data
#' @keywords datasets
#' @name ADHD
#' @usage data("ADHD")
#' @references
#' \insertAllCited{}
NULL

#' Short Boredom Proneness Scale Responses
#'
#' This dataset includes responses to the 8-item Short Boredom Proneness Scale (SBPS),
#' a self-report measure of an individual's susceptibility to boredom
#' \insertCite{Martarelli_2023_Boredom}{bgms}. Items were rated on a 7-point Likert scale
#' ranging from 1 ("strongly disagree") to 7 ("strongly agree"). The scale was administered
#' in either English \insertCite{Struk_2015_boredom}{bgms} or French (translated by \insertCite{Martarelli_2023_Boredom}{bgms}).
#'
#' @format A data frame with 986 rows and 9 columns. Each row corresponds to a respondent.
#' \describe{
#'   \item{language}{Language in which the SBPS was administered: "en" = English, "fr" = French}
#'   \item{loose_ends}{I often find myself at "loose ends," not knowing what to
#'   do.}
#'   \item{entertain}{I find it hard to entertain myself.}
#'   \item{repetitive}{Many things I have to do are repetitive and monotonous.}
#'   \item{stimulation}{It takes more stimulation to get me going than most
#'   people.}
#'   \item{motivated}{I don't feel motivated by most things that I do.}
#'   \item{keep_interest}{In most situations, it is hard for me to find
#'   something to do or see to keep me interested.}
#'   \item{sit_around}{Much of the time, I just sit around doing nothing.}
#'   \item{half_dead_dull}{Unless I am doing something exciting, even dangerous,
#'   I feel half-dead and dull.}
#' }
#'
#' @source \insertCite{Martarelli_2023_Boredom;textual}{bgms}.
#' Data retrieved from \url{https://osf.io/qhux8}.
#' Licensed under the CC-BY 4.0: \url{https://creativecommons.org/licenses/by/4.0/}
#'
#'
#' @docType data
#' @keywords datasets
#' @name Boredom
#' @usage data("Boredom")
#' @references
#' \insertAllCited{}
NULL
