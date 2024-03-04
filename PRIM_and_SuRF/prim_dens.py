import numpy as np

class PRIMdens:
    def __init__(self, X, alpha=0.05):
        self.alpha = alpha
        self.boxes_ = []  # To store intermediate boxes
        self.densities_ = []  # To store densities of intermediate boxes
        self.X_ = X  	

    def fit(self):
        # Initial box is the unit box
        box = np.vstack((np.zeros(self.X_.shape[1]), np.ones(self.X_.shape[1])))
        n_points = self.X_.shape[0]
        
        for iteration in range(120):
            max_density = -np.inf
            best_cut = None
            best_box = None
            n_in_best_box = n_points  # Initially, all points are inside the box
            ind_in_best_box = None # Initially, all points are inside the box

            for dim in range(self.X_.shape[1]):
                for direction in [0, 1]:  # 0 for lower bound, 1 for upper bound
                    trial_box = box.copy()
                    current_dim_length = trial_box[1, dim] - trial_box[0, dim]
                    
                    # Adjust the boundary in the dimension such that the resulting box's volume is (1-alpha) times the previous box
                    adjustment = current_dim_length * self.alpha
                    if direction == 0:
                        trial_box[0, dim] += adjustment
                        in_box = self.X_[:,dim] >= trial_box[0, dim]
                    else:
                        trial_box[1, dim] -= adjustment
                        in_box = self.X_[:,dim] <= trial_box[1, dim]
                    
                    # Count points in the trial box
                    n_in_box = np.count_nonzero(in_box)
                    
                    # Calculate volume of the trial box
                    volume = np.prod(trial_box[1] - trial_box[0])
                    
                    # Calculate density
                    density = n_in_box / volume if volume > 0 else 0
                    
                    if density > max_density:
                        max_density = density
                        best_cut = (dim, direction)
                        best_box = trial_box
                        n_in_best_box = n_in_box  # Update the count for best box
                        ind_in_best_box = in_box

            # If we can't find a cut that improves density or the box contains less than 1/alpha points, break
            if best_cut is None or n_in_best_box < 1/self.alpha:
                break

            # Update the box and record it
            box = best_box
            self.boxes_.append(box)
            self.densities_.append(max_density)
            self.X_ = self.X_[ind_in_best_box]

        return self

    def get_boxes(self):
        return self.boxes_

    def get_densities(self):
        return self.densities_