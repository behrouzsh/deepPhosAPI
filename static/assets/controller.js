/**
 * Created by shamsabz on 6/14/18.
 */
var appModule = angular.module('projectModule', ['angular.filter', 'angularUtils.directives.dirPagination', 'mgcrea.ngStrap']);



//function appModule($scope) {

//$scope.exportItems = [{a:1,b:10},{a:2,b:20},{a:3,b:30}];
//};

//This is for being able to download PLN
appModule.config(['$compileProvider',
    function ($compileProvider) {
        $compileProvider.aHrefSanitizationWhitelist(/^\s*(https?|ftp|mailto|tel|file|blob):/);
    }]);

appModule.factory('SharedService', function( $http, $q, $location) {
    //Defining Start point variables
    var self = this;
    var vars = {};
    vars.textArea = "IYQY[+80]IQSR[+42] 1.3\n[aK]AYSFCGTVE[pY]MAP -2.14\nKAY[+80]SF[myC]GTVE[pY]MAPEVVNR[+42.1] -0.5";


    var deferred = $q.defer();


    return {

        getSiteVisit : function() {
            return $http.get("api/increment/" + 0)
                .then(function(response) {
                    // promise is fulfilled
                    console.log(response);
                    deferred.resolve(response.data);
                    return response.data;
                    //return response.data;
                }, function(response) {
                    // the following line rejects the promise
                    deferred.reject(response);
                    return deferred.promise;
                });
        },

        getVar: function (variable) {
            return vars[variable.toString()];
        },

        setVar: function (key,value) {
            vars[key.toString()] = value;
        }
    }

});

appModule.controller("AboutCtrl", ['$scope', '$http', '$location', '$window', '$timeout', '$routeParams', '$filter', '$q', 'filterFilter', 'SharedService', function ($scope, $http, $location, $window, $timeout, $routeParams, $filter, $q, filterFilter, SharedService) {
    console.log("--------------- Restarting About! ---------------");

}]);
