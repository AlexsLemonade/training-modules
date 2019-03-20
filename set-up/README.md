### Prepare the RNA-seq Environment using Kitematic

1. Pull the appropriate image using command line.

- In *Mac*, search for and open `Terminal`.
- In *Windows*, search for and open `Command Prompt`.

  In your respective command line interface, copy and paste the following:
```
docker pull ccdl/training_rnaseq:2019-houston
```

2. Run the container. Change the `<PASSWORD>` in the line below to whatever you'd
  like.
```
docker run -e PASSWORD=<PASSWORD> -p 8787:8787 ccdl/training_rnaseq:2019-houston
```

3. Open `Kitematic` - you should see an image running.

4. `Settings` > `Volumes` > Set local folder to "main directory" that was
transferred from the flash/hard drive, using the `CHANGE` button.
![Folder](screenshots/all-02-volume.png)

5. Navigate to RStudio window.

  - In a *Windows* or *Mac* in Kitematic, go to the `Settings` > `Hostname/Ports`
    tab and click on the blue lettering.
![Folder](screenshots/all-01-network.png)

  - Alternatively, for a *Mac*, you can navigate to the RStudio window by typing
    `localhost:8787` in your web browser

6. Log into `RStudio`. The username will be 'rstudio' and the password will be
whatever you selected above (can also be accessed from the `Settings` >
`General` panel).

7. You should see a `kitematic/` folder in your `RStudio` Files panel. Click on it.
If you do not see the training modules folders in the kitematic folder, raise
your hand.

8. Go to Docker preferences.
  - In a *Mac*, click the Docker whale icon in your task bar and go to `Preferences.`
  ![Folder](screenshots/preferences_tab.png)
  - In a *Windows*, right click on Docker application and go to `Settings`.

9. Set Docker resources limits.
Go to the `Advanced` panel. This contains options that allow you to set how
much of your computer's resources Docker is allowed to use.
You probably want to leave at least 2GB of your computer's Memory that isn't
available to Docker.  
  ![Folder](screenshots/advanced_panel.png)

### If Kitematic doesn't work:

If all else fails and Kitematic is not working for you, go to your `Terminal` or
`Command Prompt` (for Mac or Windows respectively) and type in the following, but
replacing <PATH_TO_TRAINING_FOLDERS> bit with the absolute path to
"main directory" that was transferred from the flash/hard drive.
```
docker run -it --rm --mount type=volume,dst=/home/rstudio/kitematic,volume-driver=local,volume-opt=type=none,volume-opt=o=bind,volume-opt=device=<PATH_TO_TRAINING_FOLDERS> -e PASSWORD=<PASSWORD> -p 8787:8787
```
After starting your container this way, you can get to the RStudio window in
a similar way as described above:
- In Mac, type: `localhost:8787` in your web browser.
- In Windows, go to `Command Prompt`, type: `ipconfig` and click enter.
  Find the number that corresponds to the `Virtual Box Host Network` and the
  `IPv4 Address`. Copy and paste it.
  Put that number and `:8787` at the end of it in your browser.

Resume with step 6 and 7.
