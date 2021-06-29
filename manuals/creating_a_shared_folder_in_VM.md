# Setting u a shared folder between your computer and your VM

Once you start working on your own data after the course, you may want to import some of this data from your own computer to the VM to work on. This will allow you to make use of the tools you now have some experience with.

1. Create a directory to share called ‘VMshare’ (or whatever you want to call it) on your own machine. I suggest you create another directory within this first one where you will do your work.
2. Open VirtualBox Manager, but don't load your VM. Select the VM (without starting it), and click ‘Settings’ in the top menu bar. 
3. Go to ‘Shared Folders’ and select the ‘+’ button on the right to "Add a shared folder". 
4. In the ‘Folder Path’ select ‘Other’ and navigate to and select the ‘VMshare’ folder that you have created. Check the box "Auto-mount", and then click on ‘OK’. 
5. Click "ok" again, and then start your VM as normal. Once loaded, you should see a directory on your desktop showing the shared space. It may ask you for a password, which should be "manager".
6. Finally, open a terminal and run 'sudo usermod -a -G vboxsf manager'. 
7. This should complete the process and allow you t share between your computer and the VM.

Caveat
- VM has a 100GB hard disk so you're limited to this, so be conscious of what you are loading
- If you have "big data", you might consider other options for analysing your data rather than the VM, such as a local high performance computer setup or facility. The VM running on yourown computer will not work that well for this type of data.
