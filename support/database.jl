using MySQL

path = "../data/database/mysql/"

"""
CREATE TABLE `FBD8` (
	`id` INT NOT NULL AUTO_INCREMENT,
	`seed` INT NOT NULL,
	`distribution` VARCHAR(255) NOT NULL DEFAULT 'uniform',
	`loadSharing` VARCHAR(255) NOT NULL,
	`t` FLOAT NOT NULL,
	`alpha` FLOAT NOT NULL,
	`beta` FLOAT NOT NULL,
	`gamma` FLOAT NOT NULL,
	`lastStep` INT NOT NULL,
	`simulationTime` FLOAT NOT NULL,
	`breakSequence` blob(256) NOT NULL,
	`nrClusters` blob(256) NOT NULL,
	`largestCluster` blob(256) NOT NULL,
	`largestPerimeter` blob(256) NOT NULL,
	`mostStressedFiber` blob(256) NOT NULL,
	`spanningClusterSize` INT NOT NULL,
	`spanningClusterPerimeter` INT NOT NULL,
	`spanningClusterStep` INT NOT NULL,
	`spanningCenterOfMass` blob(512) NOT NULL,
	`spanningRadiusOfGyration` blob(512),
	PRIMARY KEY (`id`)
);

CREATE TABLE `FBD1024` (
	`id` INT NOT NULL AUTO_INCREMENT,
	`seed` INT NOT NULL,
	`distribution` VARCHAR(255) NOT NULL DEFAULT 'uniform',
	`loadSharing` VARCHAR(255) NOT NULL,
	`t` FLOAT NOT NULL,
	`alpha` FLOAT NOT NULL,
	`beta` FLOAT NOT NULL,
	`gamma` FLOAT NOT NULL,
	`lastStep` INT NOT NULL,
	`simulationTime` FLOAT NOT NULL,
	`breakSequence` blob(4194304) NOT NULL,
	`nrClusters` blob(4194304) NOT NULL,
	`largestCluster` blob(4194304) NOT NULL,
	`largestPerimeter` blob(4194304) NOT NULL,
	`mostStressedFiber` blob(4194304) NOT NULL,
	`spanningClusterSize` INT NOT NULL,
	`spanningClusterPerimeter` INT NOT NULL,
	`spanningClusterStep` INT NOT NULL,
	`spanningCenterOfMass` blob(8388608) NOT NULL,
	`spanningRadiusOfGyration` blob(8388608),
	PRIMARY KEY (`id`)
);

"""