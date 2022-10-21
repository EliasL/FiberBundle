using Logging
# Info is log level 0
rootLog = Logging.LogLevel(0) # Logging inside the root node
nodeLog = Logging.LogLevel(rootLog.level -1) # Logging inside a node
settingLog = Logging.LogLevel(nodeLog.level -1) # Logging inside a setting
threadLog = Logging.LogLevel(settingLog.level -1) # Logging inside a thread
seedLog = Logging.LogLevel(threadLog.level -1)
singleBreakLog = Logging.LogLevel(seedLog.level -1)
clusterLog = Logging.LogLevel(singleBreakLog.level -1)