Research Utilities
===

This repository is a collection of utility scripts and (small) programs that are used to aid in the research process.
Notable examples are:

1. pubmedAbsExtract.py
   - Utility to download abstracts from pubmed database. This can be run through a weekly cron job. The file can also be sent through e-mail via mutt
2. connectVpn.sh
   - Shell script to connect to University of Colorado VPN server. Since migrating from CISCO to globalprotect, the connection has been pretty slow... Regardless if you are using Linux then it's nice to have this start-up at the beginning of your server startup. This can also be run through a daily cron job (if need be). 
3. sendMail.sh
   - Shell script to send mail via mutt. The setup I have is just sending a .csv file that summarizes the title and abstract content that is pulled from the pubmedAbsExtract.py file above. This can also be set up as a cron job.

Contact
===

[e-mail](mailto:chong.kim@ucdenver.edu). Let me know if these scripts need updating or something else is on your mind!
