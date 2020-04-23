## RStudio Server Set Up

We will send you your Username and temporary password by a direct message on Slack, after you are logged in on [Cancer Data Science](http://ccdatalab.org/slack).

Go to https://rstudio.ccdatalab.org and type in your username and temporary password, and click `Sign in`.
You may want to bookmark https://rstudio.ccdatalab.org for the duration of the workshop.

<img src = "screenshots/rstudio-server-login.png" width = 600>

Signing in should bring you to the RStudio session page.
Click on the `Terminal` tab.

<img src = "screenshots/rstudio-session-terminal.png" width = 600>

Type in the `passwd` command in the `Terminal` tab.

<img src = "screenshots/rstudio-password-change-1.png" width = 400>

...and press `Enter`. Then type in the password you were given in Slack.

<img src = "screenshots/rstudio-password-change-2.png" width = 400>

Type in the new password you've chosen once, and press `Enter`, then it will ask you to confirm by typing it in again.
Also press `Enter`.
If you forget your password at anytime, Slack on of the CCDL team members to assist you.

## Navigating the RStudio Server

In the upper right corner of the session page, you should see these buttons:

<img src = "screenshots/rstudio-session-buttons.png" width = 400>

You can restart your R session with the orange, circular, on/off button.
You will want to do this each time you switch to notebooks, more on this in the
[RStudio guide instructions](../intro-to-R-tidyverse/00a-rstudio_guide.md).

### About RStudio Server sessions

Clicking on the house button will bring you to the RStudio server workspaces page.
Here, you will be able to see and manage all your currently running sessions.

<img src = "screenshots/rstudio-workspaces.png" width = 600>

You can start a new Session with the `+ Session` button.
We kindly ask you shut down any sessions you aren't using so we conserve computing power.
Click on a Session in the list to return to it.

Back on the session page, you are able to see the current files in your `home` directory in the `Files` tab in the lower right panel in your session.

<img src = "screenshots/rstudio-files.png" width = 600>

Here, you will find the `training-modules` folder that contains our course materials and the `shared-data` that contains the data we will be using in the modules.
The files in these folders are accessible in each R session you start.
Restarting R sessions will refresh what is in your `Environment` tab in the upper right panel.
We go into more detail on the R environment and other RStudio navigating tidbits in our [guide to RStudio](../intro-to-R-tidyverse/00a-rstudio_guide.md) as well as our [first intro to R notebook](../intro-to-R-tidyverse/01-intro_to_base_R.Rmd).

As always, please reach out to our CCDL team through Slack if you have any questions!
