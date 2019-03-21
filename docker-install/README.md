# Installing Docker

Follow the instructions below for Mac and Windows operating systems.
When you are finished, return to the schedule for your workshop for the link to the next steps.

## Pre-install Steps

1. If you do not have an account with the Docker Store, visit the [Docker Store website](https://store.docker.com/). In the upper right corner, select Log In. Create an account.

## Installation

* [Mac OS](#macos)
* [Windows 10 Pro](#windows-10-pro)
* [Windows 7](#windows-7)

### MacOS

1. Visit the [Docker Store](https://store.docker.com/search?type=edition&offering=community) and download Docker Community Edition for Mac.
![Mac Download](screenshots/mac-00-download.png)
2. Open the downloaded file and drag the Docker Whale icon to your Applications folder.
![Mac Docker Install](screenshots/mac-01-applications.png)
3. Open the Applications folder and find and open Docker.
![Mac Docker Start](screenshots/mac-02-opendocker.png)
4. Click the whale icon in the taskbar. Log in to the Docker store.
![Mac Docker Log In](screenshots/mac-03-dockerlogin.png)
5. Click the whale icon in the taskbar and go to Preferences.
![Mac preferences](screenshots/mac-preferences.png)
6. Increase the resources available to Docker. For reference, we've found that 2 CPUs and 6 GB of RAM on a quad-core with 8 GB of RAM total works for the workshop offerings. If you are able to provide more resources, do so.  
![Mac advanced](screenshots/mac-advanced-settings.png)
7. Click the whale icon in the taskbar. Click the Kitematic menu entry.  
![Mac Kitematic Round 1](screenshots/mac-04-kitematic.png)
8. Follow the prompts to download Kitematic.  
![Mac Kitematic Download](screenshots/mac-05-kitematicinstall.png)
9. Drag the Kitematic icon to the Applications folder.  
![Mac Kitematic Applications](screenshots/mac-06-kitematicapps.png)
10. Open the new Kitematic program. Accept the prompt.  
![Mac Kitematic Open](screenshots/mac-07-openkite.png)

### Windows 10 Pro

Here are the [Windows install instructions](https://docs.docker.com/docker-for-windows/install/).
We summarize the most important steps below, but if you run into trouble you may need to consult the Docker documentation.

1. Visit the [Docker Store](https://store.docker.com/search?type=edition&offering=community) and select Docker Desktop for Windows.  
![Windows 10 Download](screenshots/win10-00-download.png)
2. Download Docker Community Edition for Windows and follow the prompts.  
![Windows 10 Get Docker](screenshots/win10-01-getdocker.png)
3. During installation, a configuration menu will come up. Do not select "Use Windows containers instead of Linux containers."  
![Windows 10 configuration](screenshots/win10-configuration.png)
4. Log back in and click on the Docker Desktop icon to run Docker Desktop. Enter your Docker account credentials.    
5. Find the Docker whale in your taskbar and right click it to bring up the menu.  
![Windows Taskbar](screenshots/win10-taskbar-whale.png)
6. Go to Settings.  
![Windows Settings](screenshots/win10-taskbar-settings.png)
7. Go to the Advanced pane and increase the resources available to Docker.
For reference, we've found that 2 CPUs and 6 GB of RAM on a quad-core with 8 GB of RAM total works for the workshop offerings.  
![Windows advanced](screenshots/win10-advanced-settings.png)
7. Select Kitematic in the Docker menu.  
![Windows Taskbar Menu](screenshots/win10-taskbar-menu.png)
8. This will bring up a prompt to download Kitematic. Click download.  
![Windows download kitematic](screenshots/win10-download-kitematic.png)
9. Extract all files from the `Kitematic-Windows.zip` into `C:\Program Files\Docker\Kitematic`. You may need to create a new `Kitematic` Folder.  
10. Navigate back to the Taskbar, Docker whale, and click on Kitematic. You should now see a prompt for your Docker credentials.  

#### Troubleshooting

When you first try to run Docker Desktop, you may see an error that Hyper-V is required or that virtualization must be enabled (see [documentation](https://docs.docker.com/docker-for-windows/troubleshoot/#virtualization)).  

![virtualization error](screenshots/win10-virtualization-error.png)

You can check if you have Hyper-V support in the Performance tab of Task Manager.  

![Windows Task Manager](screenshots/win10-taskmanager-hyperv.png)

If you need to enable virtualization, you will need to enter the BIOS setup menu by restarting your computer and booting safe mode. You may find [this article](https://www.laptopmag.com/articles/access-bios-windows-10) helpful.
In the BIOS set up, you will need to enable virtualization.
Unfortunately, the steps for enabling virtualization are manufacturer- and sometimes model-specific.  

### Windows 7

1. Visit the [Docker Toolbox Docker documentation page](https://docs.docker.com/toolbox/toolbox_install_windows/). Download Docker Toolbox.  
![Windows 7 Download](screenshots/win7-00-download.png)
2. After the download completes, follow the prompts to complete the Docker Toolbox installation.  
3. After the installation, select the kitematic icon from the desktop.  
![Windows 7 kitematic](screenshots/win7-01-kitematic.png)
4. Select the option to start in VirtualBox. Then you should see a loading image while the Virtualbox session starts.  
![Windows 7 Virtual Box](screenshots/win7-02-startingvm.png)


### Now return to the workshop schedule for the link to the next steps.
