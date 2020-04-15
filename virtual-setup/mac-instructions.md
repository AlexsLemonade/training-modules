# Mac Setup Instructions for Virtual Workshops

TOC

## Zoom Installation and Setup

### New Installation
If you do not have Zoom installed yet, you will need to download the client.
To do so, go to https://zoom.us/download and click the button to download the latest version of **Zoom Client for Meetings**

![Download](screenshots/mac-zoom-01-download.png)

The `Zoom.pkg` installer package will most likely be saved to your Downloads folder.
Open the installer and follow the directions to install Zoom.
You will likely need to enter your (or an administrator) password to allow the installation to proceed.

![Zoom installer window](screenshots/mac-zoom-03-installer.png)

When the installation is completed, Zoom should automatically open, and you will see the application splash screen.
You do *not* need to login or create a Zoom user.

![Zoom splash screen](screenshots/mac-zoom-04-zoomsplash.png)

If you are on Mac OS 10.15 (Catalina), you may also be asked to give Zoom permission to access your downloads folder.
Annoyingly, this window may be blocked by the splash screen, but if you see it flash by, you will want to find it and click **OK**.

![Download folder permissions](screenshots/mac-zoom-05-download-permissions.png)

### Set Up Preferences

Open the Zoom Preferences window from the **zoom.us** menu.

![Zoom preferences menu](screenshots/mac-zoom-06-prefs-menu.png)

You may be asked right away to grant Zoom permission to use your microphone, which you should grant (click **OK**).
Again, this window may get hidden behind the Zoom application, so be on the lookout for it.

![Microphone permissions](screenshots/mac-zoom-07-mic-permissions.png)

If it is not selected already, click on the **Audio** panel in the left sidebar.
If you have granted permission for Zoom to use the microphone, you should now see green and/or red bars in the *Input Level* section as you talk.
If you want to further test your microphone and speaker settings, you can click the *Test Speaker* and *Test Mic* buttons.

We recommend that you select the checkbox to "Join audio by computer when joining a meeting" to save you the future annoyance of having to do this every time you join a meeting.

![Microphone preferences](screenshots/mac-zoom-08-mic-prefs.png)

Select **Video** from sidebar to make sure your camera is set up; you will likely have to grant permissions once again.
(Have you noticed a recurring theme?)

![Video preferences](screenshots/mac-zoom-09-video-prefs.png)

### Screen sharing permissions

MacOS has gotten stricter and stricter about security, so you will have to grant a few more permissions from the **System Preferences** application to allow for screen sharing to work.
Open **System Preferences** from the Apple menu or the Dock, and choose **Security & Privacy**, then click on the **Privacy** button at the top.
Click on the lock icon in the lower left and enter your administrator password to allow changes.
![Security preferences](screenshots/mac-zoom-10-security-prefs.png)

In the lect sidebar, select **Accessibility**, and make sure the checkbox next to **zoom.us** in the "Allow the aps below to control your computer" list is selected, as shown.

![Accessibility permissions](screenshots/mac-zoom-11-accessibility-prefs.png)

If the Zoom Application does not appear in the list, click the "**+**" button  below the list and open the application from the dialog box that appears (it should be in your *Applications* folder).

![Accessibility selection](screenshots/mac-zoom-11a-accessibility-select.png)

There is one more place where permissions will need to be given to Zoom to allow screen sharing, and that is the **Screen Recording** section of the same preferences pane.  

![Screen recording Preferences](screenshots/mac-zoom-12-recording-prefs.png)

You will probably have to wait for your first Zoom meeting to enable those permissions, as Zoom may not yet appear in that list, and there is no way to add it on your own.
Our setup meetings are a great place to test this out.
The first time you try to share your screen from Zoom, you will see a popup like the one below.
Click **OK**
![Recording popup](screenshots/mac-zoom-13-recording-popup.png)

You will then get a prompt that the Zoom app needs to be quit before this setting can be applied.
Select **Quit Now**.

![Recording restart prompt](screenshots/mac-zoom-14-recording-restart.png)

Once Zoom has quit, you can click a meeting link to reopen it, rejoining the same meeting you were in, or any other meeting.

You should now be all set to share your screen and fully participate in Zoom meetings and help sessions!



### Slack